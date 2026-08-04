package test.java.db;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import main.java.db.DBGear;
import main.java.input.CParameter;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.HashSet;

/**
 * Protein N-terminal methionine excision must not fire under the {@code NoCut} pass-through.
 *
 * <p>{@code NoCut} is how the Osprey peptide-level FASTA is predicted: every entry is ALREADY a
 * peptide, so there is no protein N-terminus to excise from. Before this guard the clip block in
 * {@link DBGear#digest_protein} fired anyway, because its "peptide is a prefix of the protein"
 * filter is trivially true when the entry IS the protein - so every M-initial peptide silently
 * gained a clipped copy that no digest ever produced.</p>
 *
 * <p>On the peptide-level FASTA that is not cosmetic. A target and its entrapment are separate
 * entries, and the entrapment shuffle preserves only the C-terminus, so M-initial status is
 * uncorrelated within a quartet and the clip fires on one side only. Measured on the 1,391,588
 * quartet Astral build: 24,093 entrapment peptides were M-initial where their target was not
 * (19,560 of them still inside the 400-900 m/z window, which reproduced the 19,559 orphans Osprey
 * reported and dropped), and 45,537 targets were M-initial where their entrapment was not, leaving
 * those targets with no entrapment coverage at all.</p>
 */
public class NoCutMetClipTest {

    private static final String M_INITIAL_PEPTIDE = "MPEPTIDEK";

    private int savedEnzyme, savedMissed, savedMinLen, savedMaxLen;
    private boolean savedClipM;

    @BeforeMethod
    public void saveCParameter() {
        DBGear.init_enzymes();
        savedEnzyme = CParameter.enzyme;
        savedMissed = CParameter.maxMissedCleavages;
        savedMinLen = CParameter.minPeptideLength;
        savedMaxLen = CParameter.maxPeptideLength;
        savedClipM = CParameter.clip_nTerm_M;

        CParameter.maxMissedCleavages = 1;
        CParameter.minPeptideLength = 7;
        CParameter.maxPeptideLength = 35;
        CParameter.clip_nTerm_M = true;
    }

    @AfterMethod
    public void restoreCParameter() {
        CParameter.enzyme = savedEnzyme;
        CParameter.maxMissedCleavages = savedMissed;
        CParameter.minPeptideLength = savedMinLen;
        CParameter.maxPeptideLength = savedMaxLen;
        CParameter.clip_nTerm_M = savedClipM;
    }

    @Test
    public void testNoCutDoesNotClipNTerminalMethionine() {
        // Looked up by literal name rather than the constant added with the fix, so this test
        // compiles against the unfixed tree and can be watched to fail there.
        Enzyme noCut = DBGear.getEnzymeByIndex(DBGear.getEnzymeIndexByName("NoCut"));
        HashSet<String> peptides = new DBGear().digest_protein(noCut, M_INITIAL_PEPTIDE);

        Assert.assertTrue(peptides.contains(M_INITIAL_PEPTIDE),
                "NoCut must pass the entry through unchanged");
        Assert.assertFalse(peptides.contains(M_INITIAL_PEPTIDE.substring(1)),
                "NoCut must NOT mint a Met-clipped copy: the entry is already a peptide, so there " +
                        "is no protein N-terminus to excise, and the clip fires asymmetrically " +
                        "across a target/entrapment quartet");
        Assert.assertEquals(peptides.size(), 1,
                "NoCut on a single peptide must yield exactly that peptide");
    }

    /**
     * The other half of the property. Without this, the guard could be satisfied by disabling the
     * clip everywhere, which would drop legitimate protein N-terminal peptide forms from every
     * ordinary digest.
     */
    @Test
    public void testRealEnzymeStillClipsNTerminalMethionine() {
        Enzyme trypsin = DBGear.getEnzymeByIndex(DBGear.getEnzymeIndexByName("Trypsin"));
        HashSet<String> peptides = new DBGear().digest_protein(trypsin, "MPEPTIDEKSAMPLERPEPTIDEK");

        Assert.assertTrue(peptides.contains("MPEPTIDEK"),
                "the unclipped protein N-terminal peptide must still be produced");
        Assert.assertTrue(peptides.contains("PEPTIDEK"),
                "a genuine protein N-terminus must still yield its Met-clipped form");
    }
}
