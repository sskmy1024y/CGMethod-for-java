import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import org.la4j.*;

public class CG {
  private static double EPSILON = 1e-10;
  private static long STARTTIME;
  private static long ENDTIME;

  public static void main(String[] args) {

    Matrix A = Matrix.fromCSV(loadCSV(args[0]));

    System.out.println("input A ----------------");
    System.out.print(A);
    System.out.println("------------------------\n");

    Vector B = Vector.fromCSV(loadCSV(args[1]));

    System.out.println("input b ----------------");
    System.out.println(B);
    System.out.println("------------------------\n");

    double[] x0 = new double[A.rows()];
    Arrays.fill(x0, 0.0);
    
    System.out.println("Solve using CG method...\n");

    System.out.println("Result x ---------------");
    System.out.println(resolveCG(Vector.fromArray(x0), A, B));
    System.out.println("------------------------\n");

    System.out.println("Execution time: " + (ENDTIME - STARTTIME) + "ms");

  }

  public static Vector resolveCG(Vector x0, Matrix A, Vector B) {
    STARTTIME = System.currentTimeMillis();

    Vector x = x0.copy();
    Vector r = B.subtract(A.multiply(x0));
    Vector p = r.copy();
    double a, b;
    Vector r0;

    Vector Ap;

    for (int i = 0; i < 100000000; i++) {
      Ap = A.multiply(p);
      r0 = r.copy();
      // ① a
      a = r.innerProduct(r) / p.innerProduct(Ap);
      // ② x[a+1]
      x = x.add(p.multiply(a));
      // ③
      r = r.subtract(Ap.multiply(a));
      // ④
      if (r.norm() < B.norm() * EPSILON)
        break;
      // ⑤
      b = r.innerProduct(r) / r0.innerProduct(r0);
      // ⑥
      p = r.add(p.multiply(b));
    }
    ENDTIME = System.currentTimeMillis();
    return x;
  }

  public static String loadCSV(String path) {
    String csv = "";
    try (BufferedReader br = new BufferedReader(new FileReader(path))) {
      String line = null;
      while ((line = br.readLine()) != null) {
        csv += line + "\n";
      }
    } catch (FileNotFoundException e) {
      System.out.println("File not found");
    } catch (IOException e) {
      e.printStackTrace();
    }
    return csv;
  }
}
