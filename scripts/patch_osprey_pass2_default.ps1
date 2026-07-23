<#
.SYNOPSIS
    Repoints Osprey's 2nd-pass q-value default to protein-compact for a Carafe
    installer build, without touching the upstream pwiz branch.

.DESCRIPTION
    Osprey selects its merge-node 2nd-pass q-value mode at process start from the
    OSPREY_PASS2_QVALUE environment variable (Osprey.Core/OspreyEnvironment.cs,
    NormalizePass2QValue). A shipped MSI is run by collaborators who never set that
    variable, so the compiled UNSET default is what actually takes effect. Upstream
    that default is "percolator" (retrain the 2nd-pass SVM). This build wants
    "protein-compact" instead: frozen 1st-pass weights + protein-stratum competition,
    no retrain.

    This script rewrites ONLY the IsNullOrWhiteSpace(raw) (unset) branch of
    NormalizePass2QValue. The separate unrecognized-token fallback is left as
    "percolator" so an explicit OSPREY_PASS2_QVALUE=percolator still resolves to
    percolator (it falls through to that fallback) and typos still warn.

    It also asserts the OSPREY_PICK_LEGACY default is still "unset -> improved LDA
    peak picking" (PickLegacy = IsSetAndNotZero(...)). The pwiz branch may run its
    golden regression with OSPREY_PICK_LEGACY=1, but that is a runtime env var in the
    test harness; nothing in the shipped binary sets it, so the MSI already uses the
    improved pick. The assertion is a tripwire: if a future branch commit hardcodes
    the legacy pick as the default, the build fails loudly here instead of silently
    shipping legacy peak picking.

    Every match is asserted to occur exactly once. Any mismatch (upstream refactor,
    already-patched source, changed default) exits non-zero so the installer build
    fails visibly rather than silently shipping the wrong behavior. The regex is
    whitespace/newline flexible so it works on either LF or CRLF checkouts.

.PARAMETER OspreyEnvPath
    Path to Osprey.Core/OspreyEnvironment.cs in the checked-out pwiz tree.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$OspreyEnvPath
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $OspreyEnvPath)) {
    Write-Error "OspreyEnvironment.cs not found at '$OspreyEnvPath'. Osprey source layout changed?"
    exit 1
}

$src = [System.IO.File]::ReadAllText($OspreyEnvPath)

# Tripwire: the shipped binary must default to the improved LDA peak pick, i.e.
# PickLegacy stays "true only when OSPREY_PICK_LEGACY is set and non-zero".
$pickLegacyPattern = 'PickLegacy\s*=\s*IsSetAndNotZero\(\s*@?"OSPREY_PICK_LEGACY"\s*\)'
if ([regex]::Matches($src, $pickLegacyPattern).Count -ne 1) {
    Write-Error ("OSPREY_PICK_LEGACY default guard failed: expected exactly one " +
        "'PickLegacy = IsSetAndNotZero(""OSPREY_PICK_LEGACY"")'. The peak-pick default may " +
        "have changed; refusing to build an installer that could ship legacy peak picking.")
    exit 1
}

# The UNSET (IsNullOrWhiteSpace) default branch of NormalizePass2QValue. Capture the
# 'if (...)' + whitespace + 'return ' prefix so its formatting/newline is preserved;
# only the returned constant is swapped. \s* spans the line break on LF or CRLF.
$unsetDefault = '(if\s*\(\s*string\.IsNullOrWhiteSpace\(\s*raw\s*\)\s*\)\s*return\s+)PASS2_QVALUE_PERCOLATOR;'
$matches = [regex]::Matches($src, $unsetDefault)
if ($matches.Count -ne 1) {
    Write-Error ("Expected exactly one unset (IsNullOrWhiteSpace) 'return " +
        "PASS2_QVALUE_PERCOLATOR;' default in NormalizePass2QValue, found $($matches.Count). " +
        "OspreyEnvironment.cs changed upstream; re-verify the 2nd-pass q-value default before building.")
    exit 1
}

$patched = [regex]::Replace($src, $unsetDefault, '${1}PASS2_QVALUE_PROTEIN_COMPACT;')

# Write UTF-8 without BOM; regex replace preserved the original newlines in $patched.
[System.IO.File]::WriteAllText($OspreyEnvPath, $patched, (New-Object System.Text.UTF8Encoding($false)))

Write-Host "Patched OSPREY_PASS2_QVALUE unset default: percolator -> protein-compact"
Write-Host "  (frozen 1st-pass weights + protein-stratum competition, no retrain)"
Write-Host "Verified OSPREY_PICK_LEGACY default unchanged (unset -> improved LDA peak picking)"
