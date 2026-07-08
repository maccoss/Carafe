package test.java.db;

import main.java.db.PairingManifestReconciler;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Coverage for {@link PairingManifestReconciler}: pruning a pairing manifest to the peptides actually
 * present in the predicted library, dropping a pair group whose target did not survive (so no
 * entrapment is left without a target), keeping a target+decoy when only the entrapment was dropped,
 * renumbering pair indices contiguously, and I&harr;L normalization.
 */
public class PairingManifestReconcilerTest {

    private Path writeTsv(String prefix, String content) throws IOException {
        Path p = Files.createTempFile(prefix, ".tsv");
        Files.writeString(p, content, StandardCharsets.UTF_8);
        return p;
    }

    private List<String[]> readManifest(Path manifest) throws IOException {
        List<String> lines = Files.readAllLines(manifest, StandardCharsets.UTF_8);
        List<String[]> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            rows.add(lines.get(i).split("\t", -1));
        }
        return rows;
    }

    @Test
    public void prunesToLibrary_dropsGroupWithoutTarget_keepsPartialGroup_renumbers() throws IOException {
        // Prelim manifest: three quartet-ish groups plus an I/L group.
        String manifest = "sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index\n"
                // group 0: everything survives
                + "TARGETONER\tNo\tsp|P1\ttarget\t0\n"
                + "ENTRAPONER\tNo\tsp|P1_p_target\tp_target\t0\n"
                + "DECOYONER\tYes\tdecoy_sp|P1\tdecoy\t0\n"
                + "PDECOYONER\tYes\tdecoy_sp|P1_p_target\tp_decoy\t0\n"
                // group 1: target + decoy survive, entrapment dropped in prediction
                + "TARGETTWOR\tNo\tsp|P2\ttarget\t1\n"
                + "ENTRAPTWOR\tNo\tsp|P2_p_target\tp_target\t1\n"
                + "DECOYTWOR\tYes\tdecoy_sp|P2\tdecoy\t1\n"
                // group 2: target dropped -> whole group must go (no orphan entrapment/decoy)
                + "TARGETTHREER\tNo\tsp|P3\ttarget\t2\n"
                + "ENTRAPTHREER\tNo\tsp|P3_p_target\tp_target\t2\n"
                + "DECOYTHREER\tYes\tdecoy_sp|P3\tdecoy\t2\n"
                // group 3: manifest uses I, library uses L -> must match under I<->L normalization
                + "PEPTIDER\tNo\tsp|P4\ttarget\t3\n"
                + "DECOYFOURK\tYes\tdecoy_sp|P4\tdecoy\t3\n";

        // Library contains only the surviving peptides (group 2's target absent; group 1's entrapment
        // absent; group 3's target spelled with L instead of I).
        String library = "ModifiedPeptide\tStrippedPeptide\tPrecursorMz\tPrecursorCharge\n"
                + "_TARGETONER_\tTARGETONER\t500\t2\n"
                + "_ENTRAPONER_\tENTRAPONER\t500\t2\n"
                + "_DECOYONER_\tDECOYONER\t500\t2\n"
                + "_PDECOYONER_\tPDECOYONER\t500\t2\n"
                + "_TARGETTWOR_\tTARGETTWOR\t500\t2\n"
                + "_DECOYTWOR_\tDECOYTWOR\t500\t2\n"
                + "_ENTRAPTHREER_\tENTRAPTHREER\t500\t2\n" // present but its target is not -> dropped anyway
                + "_PEPTLDER_\tPEPTLDER\t500\t2\n"          // L where the manifest has I
                + "_DECOYFOURK_\tDECOYFOURK\t500\t2\n";

        Path manIn = writeTsv("recon_man", manifest);
        Path lib = writeTsv("recon_lib", library);
        Path manOut = Files.createTempFile("recon_out", ".tsv");

        PairingManifestReconciler.Result r = PairingManifestReconciler.run(
                manIn.toString(), lib.toString(), manOut.toString());

        List<String[]> rows = readManifest(manOut);
        Set<String> keptSeqs = new HashSet<>();
        Set<String> keptPairIdx = new HashSet<>();
        for (String[] row : rows) {
            keptSeqs.add(row[0]);
            keptPairIdx.add(row[4]);
        }

        // Group 0 fully retained.
        Assert.assertTrue(keptSeqs.contains("TARGETONER"));
        Assert.assertTrue(keptSeqs.contains("ENTRAPONER"));
        Assert.assertTrue(keptSeqs.contains("DECOYONER"));
        Assert.assertTrue(keptSeqs.contains("PDECOYONER"));
        // Group 1: target + decoy kept, entrapment dropped (not in library).
        Assert.assertTrue(keptSeqs.contains("TARGETTWOR"));
        Assert.assertTrue(keptSeqs.contains("DECOYTWOR"));
        Assert.assertFalse(keptSeqs.contains("ENTRAPTWOR"));
        // Group 2: target absent -> entire group dropped, including a library-present entrapment.
        Assert.assertFalse(keptSeqs.contains("TARGETTHREER"));
        Assert.assertFalse(keptSeqs.contains("ENTRAPTHREER"),
                "an entrapment whose target did not survive must not be left orphaned");
        Assert.assertFalse(keptSeqs.contains("DECOYTHREER"));
        // Group 3: kept via I<->L normalization (manifest PEPTIDER == library PEPTLDER).
        Assert.assertTrue(keptSeqs.contains("PEPTIDER"), "I/L-equivalent target must be kept");
        Assert.assertTrue(keptSeqs.contains("DECOYFOURK"));

        // Pair indices renumbered contiguously from 0 (3 kept groups -> {0,1,2}).
        Assert.assertEquals(keptPairIdx, new HashSet<>(List.of("0", "1", "2")));

        // No p_target row without a target in its (new) group.
        for (String[] row : rows) {
            if (row[3].equals("p_target")) {
                boolean hasTarget = false;
                for (String[] other : rows) {
                    if (other[4].equals(row[4]) && other[3].equals("target")) {
                        hasTarget = true;
                        break;
                    }
                }
                Assert.assertTrue(hasTarget, "p_target " + row[0] + " must share a group with a target");
            }
        }

        Assert.assertEquals(r.groupsKept, 3);
        Assert.assertEquals(r.groupsDropped, 1);
        Assert.assertEquals(r.keptTargetsWithoutDecoy, 0, "every kept group here retained its decoy");
        Assert.assertEquals(r.libraryPeptidesNotInManifest, 0,
                "the test library contains only peptides present in the manifest");
    }
}
