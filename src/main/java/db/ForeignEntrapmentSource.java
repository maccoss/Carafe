package main.java.db;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

/**
 * A pool of foreign-species peptides used as entrapment sequences instead of shuffling the target.
 *
 * <p><b>Why.</b> A shuffled entrapment is an anagram of its own target: it shares the target's
 * exact residue composition and therefore many of its fragment masses, so it is over-identified
 * and the FDP it measures is too high. A real peptide from a phylogenetically distant species
 * (Arabidopsis is the standard choice) has no such relationship to the target while still being
 * a real, physically plausible peptide that is genuinely absent from a human sample. Measured on
 * Stellar, swapping shuffle for matched Arabidopsis moved measured FDP from 1.62% to 1.15%.</p>
 *
 * <p><b>The matched draw.</b> Each target is paired 1:1 with a UNIQUE foreign peptide of the
 * nearest available neutral mass, preferring the same length. Mass matching is what makes this a
 * controlled comparison: in DIA a library entry is only scored in the isolation window containing
 * its precursor m/z, so an entrapment set drawn without regard to mass would sample a different
 * difficulty regime than the targets and bias the FDP estimate. Candidates are absence-filtered
 * (no foreign peptide equal to a real target) and pass the same
 * {@link DecoySimilarityGate} as every other generated sequence.</p>
 *
 * <p><b>Determinism.</b> Draws are served from a sorted structure and consumed in a caller-fixed
 * order, so the same inputs always produce the same assignment. Nothing here uses randomness
 * except the ratio subset, which is seeded.</p>
 */
public final class ForeignEntrapmentSource {

    /**
     * Candidates examined per target before giving up on a length band. Only a candidate that
     * fails {@link DecoySimilarityGate} consumes an attempt; a foreign peptide sharing over 40%
     * of a target's ladder is rare, so this is generous.
     */
    private static final int MAX_GATE_ATTEMPTS_PER_LENGTH = 8;

    /** Length bands searched outward from the target's length before declaring no candidate. */
    private static final int MAX_LENGTH_SEARCH_BANDS = 6;

    /** length -> (neutral mass -> peptides of that mass, not yet drawn). */
    private final Map<Integer, NavigableMap<Double, ArrayDeque<String>>> byLength = new HashMap<>();
    private int available;

    private ForeignEntrapmentSource() {
    }

    /** Number of distinct foreign peptides still available to draw. */
    public int available() {
        return available;
    }

    /**
     * Digest a foreign proteome into a draw pool.
     *
     * @param fastaPath   foreign-species protein FASTA (e.g. Arabidopsis UP000006548)
     * @param enzyme      same enzyme used for the target digest
     * @param dbGear      digester (shares the global CParameter digest options)
     * @param excluded    real target sequences - a foreign peptide equal to one of these is not
     *                    absent from the sample, so it cannot serve as entrapment
     * @param applyMzFilter whether to require the peptide fit the precursor m/z window
     */
    public static ForeignEntrapmentSource build(String fastaPath, Enzyme enzyme, DBGear dbGear,
                                                Set<String> excluded, boolean applyMzFilter,
                                                int[] charges, double minMz, double maxMz)
            throws IOException {
        ForeignEntrapmentSource pool = new ForeignEntrapmentSource();
        int proteins = 0;
        int droppedHomologous = 0;
        int droppedUnknownAa = 0;
        int droppedOutOfMz = 0;
        java.util.Set<String> seen = new java.util.HashSet<>();

        File dbFile = new File(fastaPath);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        try {
            FASTAElementIterator it = reader.getIterator();
            while (it.hasNext()) {
                FASTAElement el = it.next();
                el.setLineLength(1);
                String sequence = el.getSequence().replaceAll("\\s", "").toUpperCase();
                if (sequence.isEmpty()) {
                    continue;
                }
                proteins++;
                for (String pep : dbGear.digest_protein(enzyme, sequence)) {
                    if (!seen.add(pep)) {
                        continue;
                    }
                    if (excluded.contains(pep)) {
                        droppedHomologous++;
                        continue;
                    }
                    Double mass = EntrapmentFastaGear.peptideNeutralMass(pep);
                    if (mass == null) {
                        droppedUnknownAa++;
                        continue;
                    }
                    if (applyMzFilter
                            && !EntrapmentFastaGear.fitsMzRange(mass, charges, minMz, maxMz)) {
                        droppedOutOfMz++;
                        continue;
                    }
                    pool.add(pep, mass);
                }
            }
        } finally {
            reader.close();
        }

        Cloger.getInstance().logger.info(String.format(
                "Foreign entrapment source: %d proteins -> %d unique candidate peptides "
                        + "(%d dropped equal to a real target, %d unknown AA, %d out of m/z range)",
                proteins, pool.available, droppedHomologous, droppedUnknownAa, droppedOutOfMz));
        return pool;
    }

    private void add(String pep, double mass) {
        byLength.computeIfAbsent(pep.length(), k -> new TreeMap<>())
                .computeIfAbsent(mass, k -> new ArrayDeque<>())
                .add(pep);
        available++;
    }

    /**
     * Draw the nearest-mass unused foreign peptide for {@code targetSeq}, preferring the target's
     * own length and widening the length band only when that band is exhausted. The drawn peptide
     * is removed from the pool, so every target gets a unique one. Returns {@code null} when no
     * acceptable candidate remains.
     */
    public String draw(String targetSeq, double targetMass) {
        int targetLen = targetSeq.length();
        // Length bands in order of distance: L, L-1, L+1, L-2, L+2, ...
        for (int band = 0; band < MAX_LENGTH_SEARCH_BANDS; band++) {
            for (int sign : (band == 0 ? new int[] {0} : new int[] {-1, 1})) {
                int len = targetLen + sign * band;
                if (len < 1) {
                    continue;
                }
                String hit = drawFromLength(len, targetSeq, targetMass);
                if (hit != null) {
                    return hit;
                }
            }
        }
        return null;
    }

    /**
     * Nearest-mass acceptable draw within a single length band. Candidates rejected by the
     * similarity gate are put back so a later target with a different ladder can still use them.
     */
    private String drawFromLength(int len, String targetSeq, double targetMass) {
        NavigableMap<Double, ArrayDeque<String>> band = byLength.get(len);
        if (band == null || band.isEmpty()) {
            return null;
        }
        List<String> rejected = new ArrayList<>();
        List<Double> rejectedMass = new ArrayList<>();
        String accepted = null;
        try {
            for (int attempt = 0; attempt < MAX_GATE_ATTEMPTS_PER_LENGTH; attempt++) {
                Map.Entry<Double, ArrayDeque<String>> lo = band.floorEntry(targetMass);
                Map.Entry<Double, ArrayDeque<String>> hi = band.ceilingEntry(targetMass);
                Map.Entry<Double, ArrayDeque<String>> pick;
                if (lo == null && hi == null) {
                    break;
                } else if (lo == null) {
                    pick = hi;
                } else if (hi == null) {
                    pick = lo;
                } else {
                    pick = (targetMass - lo.getKey()) <= (hi.getKey() - targetMass) ? lo : hi;
                }
                double mass = pick.getKey();
                String candidate = pick.getValue().poll();
                if (pick.getValue().isEmpty()) {
                    band.remove(mass);
                }
                available--;
                if (DecoySimilarityGate.isCandidateAcceptable(targetSeq, candidate)) {
                    accepted = candidate;
                    break;
                }
                rejected.add(candidate);
                rejectedMass.add(mass);
            }
        } finally {
            for (int i = 0; i < rejected.size(); i++) {
                add(rejected.get(i), rejectedMass.get(i));
            }
        }
        return accepted;
    }
}
