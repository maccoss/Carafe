package main.java.db;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Integrity checks on the entrapment pairing, so a broken pairing is refused at the source rather
 * than worked around downstream.
 *
 * <p>A paired FDP estimator (FDRBench, and Osprey's entrapment diagnostics) assumes every
 * entrapment peptide has a target twin. When that assumption fails the estimator either crashes or
 * silently mis-scales, and the failure surfaces hours later in a search rather than seconds after
 * library generation. Osprey grew a whole reconciliation class to absorb exactly this; these
 * checks exist so it does not have to.</p>
 *
 * <p>There are two distinct failure surfaces and they need separate checks:</p>
 * <ol>
 *   <li><b>Manifest integrity</b> ({@link #validateQuartets}) - what the generator emitted. Cheap,
 *       runs at write time, catches generator bugs.</li>
 *   <li><b>Library against manifest</b> ({@link #validateLibraryAgainstManifest}) - what the
 *       predictor produced from it. This is the one that matters most, because the library is
 *       allowed to contain peptide FORMS the manifest never had, and a form that lands on only one
 *       side of a pair is precisely the defect that is invisible to a manifest-only check.</li>
 * </ol>
 */
public final class EntrapmentPairingValidator {

    private EntrapmentPairingValidator() {
    }

    /** A single integrity violation: what was wrong, and enough examples to act on it. */
    public static final class Violation {
        public final String kind;
        public final int count;
        public final List<String> examples;

        Violation(String kind, int count, List<String> examples) {
            this.kind = kind;
            this.count = count;
            this.examples = examples;
        }

        @Override
        public String toString() {
            return kind + ": " + count + (examples.isEmpty() ? "" : " (e.g. " + String.join(", ", examples) + ")");
        }
    }

    private static final int MAX_EXAMPLES = 5;

    /**
     * Validate the quartets the generator is about to write.
     *
     * <p>Checks that every quartet has a target, that a quartet SELECTED to carry entrapment
     * actually has one (a bare target+decoy pair emitted under an entrapment design silently
     * changes the entrapment ratio), and that no sequence is claimed by two different pair
     * indices - which would make the pairing ambiguous and the FDP denominator wrong.</p>
     */
    public static List<Violation> validateQuartets(List<EntrapmentFastaGear.Quartet> quartets) {
        // Counts are tracked separately from examples: the example lists are capped, so using
        // their size as the count would under-report any violation with more than MAX_EXAMPLES
        // instances - which is exactly the case where the number matters.
        Counter missingTarget = new Counter();
        Counter selectedWithoutEntrapment = new Counter();
        Counter entrapmentEqualsTarget = new Counter();
        Counter ambiguous = new Counter();
        Map<String, Integer> seqToPair = new HashMap<>();

        int pairIdx = 0;
        for (EntrapmentFastaGear.Quartet q : quartets) {
            if (q.target == null || q.target.isEmpty()) {
                missingTarget.hit("pair_index " + pairIdx);
            }
            if (q.entrapmentSelected && q.pTarget == null) {
                selectedWithoutEntrapment.hit(q.target);
            }
            if (q.pTarget != null && q.pTarget.equals(q.target)) {
                entrapmentEqualsTarget.hit(q.target);
            }
            for (String seq : new String[] { q.target, q.pTarget, q.decoy, q.pDecoy }) {
                if (seq == null) {
                    continue;
                }
                Integer prior = seqToPair.putIfAbsent(seq, pairIdx);
                if (prior != null && prior != pairIdx) {
                    ambiguous.hit(seq);
                }
            }
            pairIdx++;
        }

        List<Violation> out = new ArrayList<>();
        missingTarget.addTo(out, "quartet with no target");
        selectedWithoutEntrapment.addTo(out, "entrapment-selected quartet with no entrapment");
        entrapmentEqualsTarget.addTo(out, "entrapment sequence equal to its target");
        ambiguous.addTo(out, "sequence claimed by two pair indices");
        return out;
    }

    /**
     * Validate a predicted library against the manifest that generated it.
     *
     * <p>The library may legitimately contain forms the manifest lacks. What it may NOT contain is
     * an entrapment peptide with no target twin: a paired estimator cannot consume that row. The
     * complementary case - a target with no entrapment twin - is not fatal, but it is a
     * subpopulation whose entrapment ratio is zero, so it silently biases any entrapment-derived
     * FDP downward and is reported too.</p>
     *
     * @param libraryEntrapment entrapment (p_target) sequences present in the predicted library
     * @param libraryTargets    target sequences present in the predicted library
     * @param manifestPairs     manifest sequence -> pair index, for every type
     * @param manifestTypes     manifest sequence -> peptide_type (target / p_target / ...)
     */
    public static List<Violation> validateLibraryAgainstManifest(
            Collection<String> libraryEntrapment,
            Collection<String> libraryTargets,
            Map<String, Integer> manifestPairs,
            Map<String, String> manifestTypes) {

        Map<Integer, String> targetOfPair = new HashMap<>();
        Map<Integer, String> entrapmentOfPair = new HashMap<>();
        for (Map.Entry<String, Integer> e : manifestPairs.entrySet()) {
            String type = manifestTypes.get(e.getKey());
            if ("target".equals(type)) {
                targetOfPair.put(e.getValue(), e.getKey());
            } else if ("p_target".equals(type)) {
                entrapmentOfPair.put(e.getValue(), e.getKey());
            }
        }
        Set<String> libTargetSet = new HashSet<>(libraryTargets);
        Set<String> libEntrapmentSet = new HashSet<>(libraryEntrapment);

        // A CLIPPED entrapment does not pair with its manifest target - it pairs with the CLIP of
        // that target, which can only exist when the target itself starts with M. Checking merely
        // that the pair has some library target is not enough, and was the first form of this
        // method: it accepted exactly the sequences that shipped broken.
        Counter orphanEntrapment = new Counter();
        for (String e : libraryEntrapment) {
            Integer idx = manifestPairs.get(e);
            String requiredTarget;
            if (idx != null) {
                requiredTarget = targetOfPair.get(idx);
            } else {
                Integer clippedIdx = manifestPairs.get("M" + e);
                if (clippedIdx == null) {
                    // Not a manifest peptide and not the Met-clip of one: an unexplained extra,
                    // equally unusable by a paired estimator.
                    orphanEntrapment.hit(e);
                    continue;
                }
                String t = targetOfPair.get(clippedIdx);
                requiredTarget = (t != null && t.startsWith("M")) ? t.substring(1) : null;
            }
            if (requiredTarget == null || !libTargetSet.contains(requiredTarget)) {
                orphanEntrapment.hit(e);
            }
        }

        Counter uncoveredTargets = new Counter();
        for (String t : libraryTargets) {
            Integer idx = manifestPairs.get(t);
            String requiredEntrapment;
            if (idx != null) {
                requiredEntrapment = entrapmentOfPair.get(idx);
            } else {
                Integer clippedIdx = manifestPairs.get("M" + t);
                if (clippedIdx == null) {
                    continue;
                }
                String e = entrapmentOfPair.get(clippedIdx);
                requiredEntrapment = (e != null && e.startsWith("M")) ? e.substring(1) : null;
            }
            if (requiredEntrapment == null || !libEntrapmentSet.contains(requiredEntrapment)) {
                uncoveredTargets.hit(t);
            }
        }

        List<Violation> out = new ArrayList<>();
        orphanEntrapment.addTo(out, "entrapment peptide with no target twin");
        uncoveredTargets.addTo(out, "target with no entrapment coverage");
        return out;
    }

    /** Render violations as a single message suitable for an exception or a log line. */
    public static String describe(List<Violation> violations) {
        StringBuilder sb = new StringBuilder();
        for (Violation v : violations) {
            if (sb.length() > 0) {
                sb.append("; ");
            }
            sb.append(v);
        }
        return sb.toString();
    }

    private static void add(List<String> list, String example) {
        if (list.size() < MAX_EXAMPLES) {
            list.add(example);
        }
    }

    /** A violation tally: full count, capped examples. */
    private static final class Counter {
        private int count;
        private final List<String> examples = new ArrayList<>();

        void hit(String example) {
            count++;
            add(examples, example);
        }

        void addTo(List<Violation> out, String kind) {
            if (count > 0) {
                out.add(new Violation(kind, count, examples));
            }
        }
    }
}
