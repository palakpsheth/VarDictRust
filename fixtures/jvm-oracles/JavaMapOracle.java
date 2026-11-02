import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Random;

/**
 * Deterministic JVM HashMap operation-log oracle generator for VarDictRust (Issues E & N).
 *
 * Generates three operation trace files executing against real OpenJDK java.util.HashMap:
 *   1. ops_put_threshold.txt — stride-16 and stride-32 collision chains for put (P),
 *      get-then-put (G), and computeIfAbsent (C) across 8 and 9 colliding keys.
 *   2. ops_getput.txt — [100,132,164], 1000..=1012, L861Q insert/remove pattern, and
 *      a 600-key sequence straddling 33554432 (1 << 25).
 *   3. ops_getput_random.txt — 80 seeded pseudo-random P/G/C/R sequences across collision
 *      strides {1, 16, 32, 64, 128, 1024, 65536}, tagged list or tree by inspecting
 *      whether any bucket contains a TreeNode.
 *
 * Usage:
 *   javac fixtures/jvm-oracles/JavaMapOracle.java
 *   java --add-opens java.base/java.util=ALL-UNNAMED -cp fixtures/jvm-oracles JavaMapOracle threshold > testdata/javamap/ops_put_threshold.txt
 *   java --add-opens java.base/java.util=ALL-UNNAMED -cp fixtures/jvm-oracles JavaMapOracle getput > testdata/javamap/ops_getput.txt
 *   java --add-opens java.base/java.util=ALL-UNNAMED -cp fixtures/jvm-oracles JavaMapOracle random > testdata/javamap/ops_getput_random.txt
 */
public class JavaMapOracle {
    private static final Field TABLE_FIELD;

    static {
        try {
            TABLE_FIELD = HashMap.class.getDeclaredField("table");
            TABLE_FIELD.setAccessible(true);
        } catch (ReflectiveOperationException e) {
            throw new ExceptionInInitializerError(e);
        }
    }

    private static int capacity(HashMap<Integer, Integer> m) {
        try {
            Object[] table = (Object[]) TABLE_FIELD.get(m);
            return table == null ? 0 : table.length;
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }

    private static boolean hasTreeNode(HashMap<Integer, Integer> m) {
        try {
            Object[] table = (Object[]) TABLE_FIELD.get(m);
            if (table == null) {
                return false;
            }
            for (Object node : table) {
                if (node != null && node.getClass().getSimpleName().contains("TreeNode")) {
                    return true;
                }
            }
            return false;
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }

    private static void emitNew(StringBuilder out, HashMap<Integer, Integer> m) {
        m.clear();
        out.append("NEW\n");
    }

    private static void emitOp(StringBuilder out, HashMap<Integer, Integer> m, char op, int k, int step) {
        switch (op) {
            case 'P':
                m.put(k, step);
                out.append("P ").append(k).append('\n');
                break;
            case 'G':
                if (m.get(k) == null) {
                    m.put(k, step);
                }
                out.append("G ").append(k).append('\n');
                break;
            case 'C':
                m.computeIfAbsent(k, x -> 0);
                out.append("C ").append(k).append('\n');
                break;
            case 'R':
                m.remove(k);
                out.append("R ").append(k).append('\n');
                break;
            default:
                throw new IllegalArgumentException("Unknown op: " + op);
        }
    }

    private static void emitDump(StringBuilder out, HashMap<Integer, Integer> m, String label) {
        int cap = capacity(m);
        boolean tree = hasTreeNode(m);
        out.append("D ").append(label)
           .append(" cap=").append(cap)
           .append(" regime=").append(tree ? "tree" : "list")
           .append(" ORDER");
        for (Integer k : m.keySet()) {
            out.append(' ').append(k);
        }
        out.append('\n');
    }

    private static void generateThreshold(StringBuilder out) {
        out.append("# generator=fixtures/jvm-oracles/JavaMapOracle.java (mode=threshold)\n");
        out.append("# java.version=").append(System.getProperty("java.version")).append('\n');

        // 1. P stride 16 (n=8 and n=9)
        HashMap<Integer, Integer> m = new HashMap<>();
        emitNew(out, m);
        int step = 0;
        for (int i = 0; i < 8; i++) {
            emitOp(out, m, 'P', i * 16, step++);
        }
        emitDump(out, m, "P_stride16_n8");
        emitOp(out, m, 'P', 8 * 16, step++);
        emitDump(out, m, "P_stride16_n9");

        // 2. P stride 32 (n=8 and n=9)
        m = new HashMap<>();
        emitNew(out, m);
        for (int i = 0; i < 8; i++) {
            emitOp(out, m, 'P', i * 32, step++);
        }
        emitDump(out, m, "P_stride32_n8");
        emitOp(out, m, 'P', 8 * 32, step++);
        emitDump(out, m, "P_stride32_n9");

        // 3. G stride 16 (n=8 and n=9)
        m = new HashMap<>();
        emitNew(out, m);
        for (int i = 0; i < 8; i++) {
            emitOp(out, m, 'G', i * 16, step++);
        }
        emitDump(out, m, "G_stride16_n8");
        emitOp(out, m, 'G', 8 * 16, step++);
        emitDump(out, m, "G_stride16_n9");

        // 4. G stride 32 (n=8 and n=9)
        m = new HashMap<>();
        emitNew(out, m);
        for (int i = 0; i < 8; i++) {
            emitOp(out, m, 'G', i * 32, step++);
        }
        emitDump(out, m, "G_stride32_n8");
        emitOp(out, m, 'G', 8 * 32, step++);
        emitDump(out, m, "G_stride32_n9");

        // 5. C stride 16 (n=8)
        m = new HashMap<>();
        emitNew(out, m);
        for (int i = 0; i < 8; i++) {
            emitOp(out, m, 'C', i * 16, step++);
        }
        emitDump(out, m, "C_stride16_n8");
    }

    private static void generateGetPut(StringBuilder out) {
        out.append("# generator=fixtures/jvm-oracles/JavaMapOracle.java (mode=getput)\n");
        out.append("# java.version=").append(System.getProperty("java.version")).append('\n');

        int step = 0;
        // 1. C vs G on [100, 132, 164]
        HashMap<Integer, Integer> m = new HashMap<>();
        emitNew(out, m);
        for (int k : new int[]{100, 132, 164}) {
            emitOp(out, m, 'C', k, step++);
        }
        emitDump(out, m, "C_100_132_164_n3");

        m = new HashMap<>();
        emitNew(out, m);
        for (int k : new int[]{100, 132, 164}) {
            emitOp(out, m, 'G', k, step++);
        }
        emitDump(out, m, "G_100_132_164_n3");

        // 2. C vs G on 1000..=1012 (13 dense keys crossing threshold 12)
        m = new HashMap<>();
        emitNew(out, m);
        for (int k = 1000; k <= 1012; k++) {
            emitOp(out, m, 'C', k, step++);
        }
        emitDump(out, m, "C_1000to1012_n13");

        m = new HashMap<>();
        emitNew(out, m);
        for (int k = 1000; k <= 1012; k++) {
            emitOp(out, m, 'G', k, step++);
        }
        emitDump(out, m, "G_1000to1012_n13");

        // 3. L861Q collision pair [55259524, 55259572] (both spread & 15 == 15)
        m = new HashMap<>();
        emitNew(out, m);
        emitOp(out, m, 'C', 55259524, step++);
        emitOp(out, m, 'C', 55259572, step++);
        emitDump(out, m, "C_L861Q_pair");

        m = new HashMap<>();
        emitNew(out, m);
        emitOp(out, m, 'G', 55259524, step++);
        emitOp(out, m, 'G', 55259572, step++);
        emitDump(out, m, "G_L861Q_pair");

        // Full L861Q window 55259401..=55259600 insert then remove all except 55259524 and 55259572
        m = new HashMap<>();
        emitNew(out, m);
        for (int k = 55259401; k <= 55259600; k++) {
            emitOp(out, m, 'G', k, step++);
        }
        for (int k = 55259401; k <= 55259600; k++) {
            if (k != 55259524 && k != 55259572) {
                emitOp(out, m, 'R', k, step++);
            }
        }
        emitDump(out, m, "G_L861Q_window");

        // 4. 600-key sequence straddling 33554432 (1 << 25)
        m = new HashMap<>();
        emitNew(out, m);
        int center = 1 << 25;
        for (int k = center - 300; k < center + 300; k++) {
            emitOp(out, m, 'G', k, step++);
        }
        emitDump(out, m, "G_straddle_33554432_n600");
    }

    private static void generateRandom(StringBuilder out) {
        out.append("# generator=fixtures/jvm-oracles/JavaMapOracle.java (mode=random)\n");
        out.append("# java.version=").append(System.getProperty("java.version")).append('\n');

        Random rng = new Random(20261101L);
        int[] strides = {1, 16, 32, 64, 128, 1024, 65536};
        int[] bases = {0, 1000, 55259400, (1 << 24) - 300, 157482700};
        int step = 0;

        for (int seq = 0; seq < 80; seq++) {
            HashMap<Integer, Integer> m = new HashMap<>();
            emitNew(out, m);
            int stride = strides[seq % strides.length];
            int base = bases[seq % bases.length];
            int span = 24 + (seq % 40);
            int nOps = 30 + (seq % 50);
            for (int i = 0; i < nOps; i++) {
                int k = base + stride * rng.nextInt(span);
                int roll = rng.nextInt(10);
                if (roll < 5) {
                    emitOp(out, m, 'G', k, step++);
                } else if (roll < 7) {
                    emitOp(out, m, 'P', k, step++);
                } else if (roll < 9) {
                    emitOp(out, m, 'C', k, step++);
                } else {
                    emitOp(out, m, 'R', k, step++);
                }
            }
            emitDump(out, m, "rand_seq_" + seq);
        }
    }

    public static void main(String[] args) {
        String mode = args.length > 0 ? args[0] : "threshold";
        StringBuilder out = new StringBuilder();
        switch (mode) {
            case "threshold":
                generateThreshold(out);
                break;
            case "getput":
                generateGetPut(out);
                break;
            case "random":
                generateRandom(out);
                break;
            default:
                throw new IllegalArgumentException("Unknown mode: " + mode);
        }
        System.out.print(out);
    }
}
