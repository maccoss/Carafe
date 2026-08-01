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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * A pool of foreign-species peptides used as entrapment sequences instead of shuffling the target.
 *
 * <p><b>Why.</b> A shuffled entrapment is an anagram of its own target: it shares the target's
 * exact residue composition and therefore many of its fragment masses, so it is over-identified
 * and the FDP it measures is too high. A real peptide from a phylogenetically distant species
 * (Arabidopsis is the standard choice) has no such relationship to the target while still being a
 * real, physically plausible peptide that is genuinely absent from a human sample.</p>
 *
 * <p><b>The matched assignment.</b> Each target is paired 1:1 with a UNIQUE foreign peptide of
 * near-identical neutral mass. Mass matching is what makes this a controlled comparison: in DIA a
 * library entry is only scored in the isolation window containing its precursor m/z, so an
 * entrapment set drawn without regard to mass would sample a different difficulty regime than the
 * targets and bias the FDP estimate.</p>
 *
 * <p><b>What the assignment optimizes, and why that choice matters.</b> The objective is to
 * MAXIMIZE THE NUMBER of pairs that co-locate - a threshold - not to minimize total mass
 * displacement. Those are different problems and they want different algorithms, which is easy to
 * get wrong:</p>
 * <ul>
 *   <li>Nearest-available in mass order accumulates a deficit (each target eats supply just above
 *       it) and every later target is dragged further off. Measured 81% co-location.</li>
 *   <li>A quantile map (rank i of n to rank i*m/n of m) is the optimal monotone transport and
 *       gives an excellent worst case, but it spreads its error evenly across every pair, so most
 *       pairs sit a few Da off. Optimal for total displacement, wrong for a threshold. Measured
 *       49%.</li>
 *   <li>Nearest-available in an order uncorrelated with mass gives most pairs a near-exact match
 *       and strands the few it serves last. Measured 95%.</li>
 * </ul>
 * <p>This implementation bins the pool finely and serves each target from the nearest non-empty
 * bin within the co-location window. Inside the window every candidate satisfies the objective
 * equally, so bins are consumed in O(1) with no search, and the preference for nearer bins keeps
 * the mass match tight as a secondary benefit. A target whose window is exhausted falls through to
 * an unbounded search, so it still gets an entrapment peptide - just not a co-locating one, which
 * is counted and reported rather than hidden.</p>
 *
 * <p><b>Determinism.</b> Bins are filled in mass-sorted order and consumed in a caller-fixed
 * order, so the same inputs always produce the same pairing.</p>
 */
public final class ForeignEntrapmentSource {

    /**
     * Neutral-mass window, in daltons, inside which an entrapment peptide is considered to
     * co-locate with its target. A 3 m/z DIA isolation window spans 6 Da at charge 2 and 9 Da at
     * charge 3; the tighter is used so a pair counted as co-locating does so at either charge.
     */
    public static final double CO_LOCATION_WINDOW_DA = 6.0;

    /** Bin width (Da). Fine enough that the nearest non-empty bin is a tight mass match. */
    private static final double BIN_WIDTH_DA = 0.25;

    /** Bins either side of the target's own bin that still count as co-locating. */
    private static final int CO_LOCATION_BINS = (int) ((CO_LOCATION_WINDOW_DA / 2.0) / BIN_WIDTH_DA);

    /** Candidates examined per target before giving up, when the gate keeps rejecting. */
    private static final int MAX_GATE_ATTEMPTS = 16;

    private final ArrayDeque<String>[] bins;
    private final double minMass;
    private int available;
    private int total;

    @SuppressWarnings("unchecked")
    private ForeignEntrapmentSource(int nBins, double minMass) {
        this.bins = new ArrayDeque[nBins];
        this.minMass = minMass;
    }

    /** Number of distinct foreign peptides still available. */
    public int available() {
        return available;
    }

    /** Pool size before any assignment. */
    public int size() {
        return total;
    }

    /**
     * Digest a foreign proteome into an assignment pool.
     *
     * @param fastaPath     foreign-species protein FASTA (e.g. Arabidopsis UP000006548)
     * @param enzyme        same enzyme used for the target digest
     * @param dbGear        digester (shares the global CParameter digest options)
     * @param excluded      real target sequences - a foreign peptide equal to one of these is not
     *                      absent from the sample, so it cannot serve as entrapment
     * @param applyMzFilter whether to require the peptide fit the precursor m/z window
     */
    public static ForeignEntrapmentSource build(String fastaPath, Enzyme enzyme, DBGear dbGear,
                                                Set<String> excluded, boolean applyMzFilter,
                                                int[] charges, double minMz, double maxMz)
            throws IOException {
        int proteins = 0;
        int droppedHomologous = 0;
        int droppedUnknownAa = 0;
        int droppedOutOfMz = 0;
        Set<String> seen = new HashSet<>();
        List<String> keptSeqs = new ArrayList<>();
        List<Double> keptMasses = new ArrayList<>();
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;

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
                    keptSeqs.add(pep);
                    keptMasses.add(mass);
                    lo = Math.min(lo, mass);
                    hi = Math.max(hi, mass);
                }
            }
        } finally {
            reader.close();
        }

        if (keptSeqs.isEmpty()) {
            throw new IOException("Foreign entrapment FASTA yielded no usable peptides: " + fastaPath);
        }
        int nBins = (int) ((hi - lo) / BIN_WIDTH_DA) + 2;
        ForeignEntrapmentSource pool = new ForeignEntrapmentSource(nBins, lo);
        // Fill in mass-sorted order so each bin's contents are deterministic and, within a bin,
        // ordered - which keeps a draw's mass error at the low end of the bin where possible.
        Integer[] order = new Integer[keptSeqs.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }
        java.util.Arrays.sort(order, (a, b) -> {
            int c = Double.compare(keptMasses.get(a), keptMasses.get(b));
            return c != 0 ? c : keptSeqs.get(a).compareTo(keptSeqs.get(b));
        });
        for (Integer idx : order) {
            pool.add(keptSeqs.get(idx), keptMasses.get(idx));
        }

        Cloger.getInstance().logger.info(String.format(
                "Foreign entrapment source: %d proteins -> %d unique candidate peptides "
                        + "(%d dropped equal to a real target, %d unknown AA, %d out of m/z range)",
                proteins, pool.total, droppedHomologous, droppedUnknownAa, droppedOutOfMz));
        return pool;
    }

    private void add(String pep, double mass) {
        int b = binOf(mass);
        if (bins[b] == null) {
            bins[b] = new ArrayDeque<>();
        }
        bins[b].add(pep);
        available++;
        total++;
    }

    private int binOf(double mass) {
        int b = (int) ((mass - minMass) / BIN_WIDTH_DA);
        if (b < 0) {
            return 0;
        }
        return Math.min(b, bins.length - 1);
    }

    /**
     * Assign a unique foreign peptide to {@code targetSeq}, preferring one that co-locates with it.
     *
     * @return the assigned peptide, or {@code null} if the pool is exhausted
     */
    public String assign(String targetSeq, double targetMass) {
        String hit = search(targetSeq, targetMass, CO_LOCATION_BINS);
        if (hit != null) {
            return hit;
        }
        // Window exhausted. Take the nearest available anywhere rather than leaving the target
        // without entrapment - the shortfall is reported as a co-location rate, not silently.
        return search(targetSeq, targetMass, bins.length);
    }

    /** Nearest non-empty bin within {@code maxBins}, serving one peptide that passes the gate. */
    private String search(String targetSeq, double targetMass, int maxBins) {
        int home = binOf(targetMass);
        int rejected = 0;
        List<String> putBack = null;
        try {
            for (int radius = 0; radius <= maxBins; radius++) {
                for (int sign = 0; sign < (radius == 0 ? 1 : 2); sign++) {
                    int b = home + (sign == 0 ? radius : -radius);
                    if (b < 0 || b >= bins.length || bins[b] == null || bins[b].isEmpty()) {
                        continue;
                    }
                    while (!bins[b].isEmpty() && rejected < MAX_GATE_ATTEMPTS) {
                        String candidate = bins[b].poll();
                        available--;
                        if (DecoySimilarityGate.isCandidateAcceptable(targetSeq, candidate)) {
                            return candidate;
                        }
                        if (putBack == null) {
                            putBack = new ArrayList<>();
                        }
                        putBack.add(candidate);
                        rejected++;
                    }
                    if (rejected >= MAX_GATE_ATTEMPTS) {
                        return null;
                    }
                }
            }
        } finally {
            // Rejected candidates go back: a different target with a different ladder can use them.
            if (putBack != null) {
                for (String s : putBack) {
                    Double m = EntrapmentFastaGear.peptideNeutralMass(s);
                    if (m != null) {
                        int b = binOf(m);
                        if (bins[b] == null) {
                            bins[b] = new ArrayDeque<>();
                        }
                        bins[b].add(s);
                        available++;
                    }
                }
            }
        }
        return null;
    }
}
