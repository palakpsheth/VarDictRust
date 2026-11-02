import com.astrazeneca.vardict.data.fishertest.FisherExact;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.Random;
import java.util.Set;

/**
 * Deterministic JVM FisherExact oracle generator for VarDictRust (Issues M & N).
 *
 * Compiles and executes against VarDict-1.8.3.jar and commons-math3-3.6.1.jar:
 *   javac -cp "$VARDICT_JAR:$CM3_JAR" fixtures/jvm-oracles/FisherOracle.java
 *   java -cp "$VARDICT_JAR:$CM3_JAR:fixtures/jvm-oracles" FisherOracle > testdata/fisher/java_fisher_oracle.tsv
 */
public class FisherOracle {
    private static void addTable(Set<String> seen, int rf, int rr, int af, int ar) {
        seen.add(rf + "\t" + rr + "\t" + af + "\t" + ar);
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.out.println("# generator=fixtures/jvm-oracles/FisherOracle.java");
        System.out.println("# vardict_java_commit=08862618f22fa54058902ac8e75ad2905b939bca");
        System.out.println("# commons_math3=3.6.1");
        System.out.println("# java.version=" + System.getProperty("java.version"));
        System.out.println("# os.arch=" + System.getProperty("os.arch"));
        System.out.println("# columns: ref_fwd\tref_rev\talt_fwd\talt_rev\tp_two_str\tp_two_bits_hex\tp_less_str\tp_less_bits_hex\tp_greater_str\tp_greater_bits_hex\todd_ratio");

        Set<String> tables = new LinkedHashSet<>();

        // 1. Explicit pin and canary tables from Issues B, G, M
        int[][] pins = {
            {0, 63, 1, 0},
            {86, 41, 0, 1},
            {35, 28, 0, 1},
            {13, 50, 0, 1},
            {145, 494, 1, 0},
            {45, 82, 1, 0},
            {60, 259, 1, 0},
            {6, 313, 1, 0},
            {9, 6, 5, 1},
            {27, 17, 3, 8},
            {10, 8, 1, 3},
            {2, 150, 7, 183},
            {14, 35, 6, 5},
            {18, 2, 12, 5},
            {25, 3, 15, 6},
            {12, 40, 8, 9},
            {0, 0, 0, 0},
            {10, 0, 0, 5},
            {0, 10, 5, 0},
            {41, 86, 1, 1},
            {7, 11, 1, 5},
        };
        for (int[] t : pins) {
            addTable(tables, t[0], t[1], t[2], t[3]);
        }

        // 2. Single-alt dyadic tie tables (N in {64, 128, 192, 256, 320, 512, 640})
        int[] dyadicTotals = {64, 128, 192, 256, 320, 512, 640};
        for (int total : dyadicTotals) {
            int refTotal = total - 1;
            int[] splits = {0, 1, 6, 13, 28, 35, 41, 45, 50, 60, refTotal / 3, refTotal / 2, refTotal - 1, refTotal};
            for (int rf : splits) {
                if (rf >= 0 && rf <= refTotal) {
                    int rr = refTotal - rf;
                    addTable(tables, rf, rr, 1, 0);
                    addTable(tables, rf, rr, 0, 1);
                }
            }
        }

        // 3. Deterministic 2x2 grid + seeded random tables to reach 250+ tables
        Random rng = new Random(424242L);
        while (tables.size() < 260) {
            int rf = rng.nextInt(120);
            int rr = rng.nextInt(120);
            int af = rng.nextInt(40);
            int ar = rng.nextInt(40);
            addTable(tables, rf, rr, af, ar);
        }

        for (String key : tables) {
            String[] parts = key.split("\t");
            int rf = Integer.parseInt(parts[0]);
            int rr = Integer.parseInt(parts[1]);
            int af = Integer.parseInt(parts[2]);
            int ar = Integer.parseInt(parts[3]);
            FisherExact fe = new FisherExact(rf, rr, af, ar);
            double pTwo = fe.getPValue();
            double pLess = fe.getPValueLess();
            double pGreater = fe.getPValueGreater();
            String orStr = fe.getOddRatio();
            System.out.printf(
                "%d\t%d\t%d\t%d\t%s\t%016x\t%s\t%016x\t%s\t%016x\t%s%n",
                rf, rr, af, ar,
                Double.toString(pTwo), Double.doubleToRawLongBits(pTwo),
                Double.toString(pLess), Double.doubleToRawLongBits(pLess),
                Double.toString(pGreater), Double.doubleToRawLongBits(pGreater),
                orStr
            );
        }
    }
}
