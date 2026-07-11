package test.java.gui;

import main.java.gui.CarafeLibraryStaging;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;

/**
 * Coverage for {@link CarafeLibraryStaging}: the local-staging publisher that keeps the Osprey
 * library steps' heavy prediction parquet (and SQLite .blib) churn off the network share, copying
 * only the finished library + report/model artifacts to the final directory.
 */
public class CarafeLibraryStagingTest {

    private Path stage;
    private Path finalDir;

    @BeforeMethod
    public void setUp() throws IOException {
        stage = Files.createTempDirectory("carafe_stage_test");
        // A sibling temp dir that does NOT yet exist, to prove publish() creates it.
        finalDir = stage.resolveSibling(stage.getFileName() + "_final");
    }

    @AfterMethod
    public void tearDown() throws IOException {
        for (Path root : new Path[] {stage, finalDir}) {
            if (root != null && Files.exists(root)) {
                try (var walk = Files.walk(root)) {
                    walk.sorted(Comparator.reverseOrder()).forEach(p -> {
                        try {
                            Files.deleteIfExists(p);
                        } catch (IOException ignored) {
                        }
                    });
                }
            }
        }
    }

    private void write(Path dir, String name) throws IOException {
        Files.write(dir.resolve(name), ("content of " + name).getBytes(StandardCharsets.UTF_8));
    }

    @Test
    public void publishesDeliverablesCreatesFinalDirAndKeepsParquetScratch() throws IOException {
        write(stage, "carafe_spectral_library.tsv");
        write(stage, "carafe_spectral_library.blib");
        write(stage, "training_report.html");
        write(stage, "rt_model.pt");
        write(stage, "peptide_forms_1.parquet");
        write(stage, "peptide_forms_11.parquet");

        int published = CarafeLibraryStaging.publish(stage.toFile(), finalDir.toFile());

        Assert.assertEquals(published, 4, "the 4 deliverables are published; the 2 parquet chunks are not");
        Assert.assertTrue(Files.exists(finalDir), "publish creates the final directory when missing");
        // Deliverables landed on the "share" (and were moved off the stage).
        for (String name : new String[] {"carafe_spectral_library.tsv", "carafe_spectral_library.blib",
                "training_report.html", "rt_model.pt"}) {
            Assert.assertTrue(Files.exists(finalDir.resolve(name)), name + " published to final dir");
            Assert.assertFalse(Files.exists(stage.resolve(name)), name + " moved out of staging");
        }
        // Prediction scratch stays on local disk and never reaches the final (network) dir.
        for (String name : new String[] {"peptide_forms_1.parquet", "peptide_forms_11.parquet"}) {
            Assert.assertTrue(Files.exists(stage.resolve(name)), name + " kept in staging");
            Assert.assertFalse(Files.exists(finalDir.resolve(name)), name + " never published");
        }
    }

    @Test
    public void overwritesAnExistingPublishedFile() throws IOException {
        Files.createDirectories(finalDir);
        Files.write(finalDir.resolve("carafe_spectral_library.tsv"),
                "stale".getBytes(StandardCharsets.UTF_8));
        write(stage, "carafe_spectral_library.tsv");

        int published = CarafeLibraryStaging.publish(stage.toFile(), finalDir.toFile());

        Assert.assertEquals(published, 1);
        Assert.assertEquals(Files.readString(finalDir.resolve("carafe_spectral_library.tsv")),
                "content of carafe_spectral_library.tsv", "REPLACE_EXISTING overwrites the stale copy");
    }

    @Test
    public void subdirectoriesAreIgnored() throws IOException {
        write(stage, "carafe_spectral_library.tsv");
        Files.createDirectories(stage.resolve("nested"));
        write(stage.resolve("nested"), "inner.txt");

        int published = CarafeLibraryStaging.publish(stage.toFile(), finalDir.toFile());

        Assert.assertEquals(published, 1, "only the top-level file is published");
        Assert.assertFalse(Files.exists(finalDir.resolve("nested")), "sub-directories are not published");
    }

    @Test(expectedExceptions = IOException.class)
    public void emptyStageThrows() throws IOException {
        CarafeLibraryStaging.publish(stage.toFile(), finalDir.toFile());
    }

    @Test(expectedExceptions = IOException.class)
    public void missingStageThrows() throws IOException {
        CarafeLibraryStaging.publish(new File(stage.toFile(), "does_not_exist"), finalDir.toFile());
    }

    @Test
    public void networkOutputPathDetectsUncSharesOnly() {
        Assert.assertTrue(CarafeLibraryStaging.isNetworkOutputPath("\\\\maccoss-nas\\home\\run"));
        Assert.assertTrue(CarafeLibraryStaging.isNetworkOutputPath("//maccoss-nas/home/run"));
        Assert.assertTrue(CarafeLibraryStaging.isNetworkOutputPath("  \\\\server\\share  "),
                "surrounding whitespace is trimmed before the UNC check");
        Assert.assertFalse(CarafeLibraryStaging.isNetworkOutputPath("C:\\Users\\maccoss\\work"));
        Assert.assertFalse(CarafeLibraryStaging.isNetworkOutputPath("D:/data/run"),
                "a local drive is not staged");
        Assert.assertFalse(CarafeLibraryStaging.isNetworkOutputPath(null));
        Assert.assertFalse(CarafeLibraryStaging.isNetworkOutputPath(""));
    }

    @Test
    public void predictionScratchPredicateMatchesOnlyPeptideFormsParquet() {
        Assert.assertTrue(CarafeLibraryStaging.isPredictionScratch("peptide_forms_1.parquet"));
        Assert.assertTrue(CarafeLibraryStaging.isPredictionScratch("peptide_forms_12.parquet"));
        Assert.assertFalse(CarafeLibraryStaging.isPredictionScratch("carafe_spectral_library.tsv"));
        Assert.assertFalse(CarafeLibraryStaging.isPredictionScratch("carafe_spectral_library.blib"));
        Assert.assertFalse(CarafeLibraryStaging.isPredictionScratch("peptide_forms_1.tsv"),
                "a non-parquet peptide_forms file is a real output, not scratch");
        Assert.assertFalse(CarafeLibraryStaging.isPredictionScratch("other.parquet"),
                "an unrelated parquet is not the prediction scratch");
    }
}
