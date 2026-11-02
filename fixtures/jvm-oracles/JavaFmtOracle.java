import com.astrazeneca.vardict.Utils;
import java.text.DecimalFormat;
import java.util.Locale;

/**
 * Deterministic JVM formatting oracle generator for VarDictRust (Issues G & M).
 *
 * Compiles and executes against VarDict-1.8.3.jar (commit 0886261) with US locale:
 *   javac -cp "$VARDICT_JAR" fixtures/jvm-oracles/JavaFmtOracle.java
 *   java -cp "$VARDICT_JAR:fixtures/jvm-oracles" JavaFmtOracle > testdata/format/java_fmt_oracle.tsv
 */
public class JavaFmtOracle {
    private static String patternForDigits(int digits) {
        if (digits <= 0) {
            return "0";
        }
        StringBuilder sb = new StringBuilder("0.");
        for (int i = 0; i < digits; i++) {
            sb.append('0');
        }
        return sb.toString();
    }

    private static void emitRow(int digits, double value) {
        String pattern = patternForDigits(digits);
        double rhe = Utils.roundHalfEven(pattern, value);
        String gvtp = Utils.getRoundedValueToPrint(pattern, value);
        String df = new DecimalFormat(pattern).format(value);
        String dts = Double.toString(value);
        long vBits = Double.doubleToRawLongBits(value);
        long rheBits = Double.doubleToRawLongBits(rhe);
        System.out.printf(
            "%d\t%s\t%s\t%016x\t%s\t%016x\t%s\t%s\t%s%n",
            digits,
            pattern,
            dts,
            vBits,
            Double.toString(rhe),
            rheBits,
            gvtp,
            df,
            dts
        );
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        System.out.println("# generator=fixtures/jvm-oracles/JavaFmtOracle.java");
        System.out.println("# vardict_java_commit=08862618f22fa54058902ac8e75ad2905b939bca");
        System.out.println("# java.version=" + System.getProperty("java.version"));
        System.out.println("# os.arch=" + System.getProperty("os.arch"));
        System.out.println("# columns: digits\tpattern\tvalue_str\tvalue_bits_hex\tround_half_even_str\tround_half_even_bits_hex\tget_rounded_value_to_print\tdecimal_format\tdouble_to_string");

        // 1. Explicit pin cases from Issue G & Issue M (D1 DecimalFormat underflow ties, trailing dots, Double.toString thresholds)
        int[] pinDigits = {
            1, 1, 1, 1, 1,
            2, 2,
            3, 3, 3, 3,
            4, 4, 4, 4, 4, 4,
            5, 5, 5,
            7, 7,
            8, 8,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4
        };
        double[] pinValues = {
            0.023, -0.023, 0.05, 2.96, 27.85,
            0.005, -0.005,
            0.0005, -0.0005, 0.0015, 1.0 / 2000.0,
            2.00001, 5.0e-5, -5.0e-5, 1.5e-4, 1.0 / 20000.0, 0.99999,
            5.0e-6, -5.0e-6, 0.015625,
            -3.45678915, -0.00000005,
            5.0e-9, -5.0e-9,
            0.0, -0.0, 1.0, 6.0, 16384.0, 9999999.0, 10000000.0, 1.0e-3, 1.0e-4, 1.2345e-4
        };
        for (int i = 0; i < pinDigits.length; i++) {
            emitRow(pinDigits[i], pinValues[i]);
        }

        // 2. Sweep +/- 5 * 10^-(d+1) across all d in 0..8
        for (int d = 0; d <= 8; d++) {
            double tie = 5.0 * Math.pow(10.0, -(d + 1));
            emitRow(d, tie);
            emitRow(d, -tie);
            emitRow(d, Math.nextUp(tie));
            emitRow(d, Math.nextDown(tie));
        }

        // 3. Double.toString scientific-notation and integer-significand boundary sweep
        double[] dtsBoundaries = {
            0.001, Math.nextDown(0.001), 0.000999999999, 0.0001, 0.0002, 0.0005,
            0.00001, 0.00005, 0.000001, 1.0e-7, 1.0e-8,
            2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 100.0, 1024.0, 16384.0, 65536.0,
            9999999.0, Math.nextDown(1.0e7), 1.0e7, 1.25e7, 1.0e8, 1.0e9
        };
        for (double v : dtsBoundaries) {
            emitRow(4, v);
            emitRow(4, -v);
        }

        // 4. Ratio grid across denominators 8, 16, 20, 40, 64, 100, 128, 200, 1000, 2000, 20000
        int[] denoms = {8, 16, 20, 40, 64, 100, 128, 200, 1000, 2000, 20000};
        int[] digitChoices = {1, 3, 4, 5};
        for (int b : denoms) {
            for (int a = 1; a <= 7; a++) {
                double ratio = (double) a / (double) b;
                for (int d : digitChoices) {
                    emitRow(d, ratio);
                }
            }
        }
    }
}
