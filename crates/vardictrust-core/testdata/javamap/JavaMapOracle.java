import java.util.*;
public class JavaMapOracle {
  public static void main(String[] a) {
    long seed = Long.parseLong(a[0]); int nOps = Integer.parseInt(a[1]);
    int base = Integer.parseInt(a[2]); int span = Integer.parseInt(a[3]); int stride = Integer.parseInt(a[4]);
    Random r = new Random(seed); HashMap<Integer, Integer> m = new HashMap<>(); StringBuilder o = new StringBuilder();
    for (int i = 0; i < nOps; i++) {
      int k = base + stride * r.nextInt(span); int op = r.nextInt(10);
      if (op < 5) { m.put(k, i); o.append("P ").append(k).append('\n'); }
      else if (op < 8) { m.computeIfAbsent(k, x -> 0); o.append("C ").append(k).append('\n'); }
      else { m.remove(k); o.append("R ").append(k).append('\n'); }
    }
    o.append("ORDER"); for (Integer k : m.keySet()) o.append(' ').append(k);
    System.out.println(o);
  }
}
