import java.util.*;

/** Quaternions. Basic operations. */
public class Quaternion {

   private double r; // real part
   private double i; // imaginary part i
   private double j; // imaginary part j
   private double k; // imaginary part k

   private static final double EPS = 1e-9;


   /** Constructor from four double values.
    * @param a real part
    * @param b imaginary part i
    * @param c imaginary part j
    * @param d imaginary part k
    */
   public Quaternion (double a, double b, double c, double d) {
      r = a;
      i = b;
      j = c;
      k = d;
   }

   /** Real part of the quaternion.
    * @return real part
    */
   public double getRpart() {
      return r;
   }

   /** Imaginary part i of the quaternion.
    * @return imaginary part i
    */
   public double getIpart() {
      return i;
   }

   /** Imaginary part j of the quaternion.
    * @return imaginary part j
    */
   public double getJpart() {
      return j;
   }

   /** Imaginary part k of the quaternion.
    * @return imaginary part k
    */
   public double getKpart() {
      return k;
   }

   /** Conversion of the quaternion to the string.
    * @return a string form of this quaternion:
    * "a+bi+cj+dk"
    * (without any brackets)
    */
   @Override
   public String toString() {
      return String.format("%.1f+%.1fi+%.1fj+%.1fk", r, i, j, k);
   }

   /** Conversion from the string to the quaternion.
    * Reverse to <code>toString</code> method.
    * @throws IllegalArgumentException if string s does not represent
    *     a quaternion (defined by the <code>toString</code> method)
    * @param s string of form produced by the <code>toString</code> method
    * @return a quaternion represented by string s
    */
   public static Quaternion valueOf (String s) {
      String[] parts = s.split("\\+");
      if (parts.length != 4) {
         throw new IllegalArgumentException("Invalid input string: " + s);
      }

      try {
         double a = Double.parseDouble(parts[0]);
         double b = Double.parseDouble(parts[1].replace("i", ""));
         double c = Double.parseDouble(parts[2].replace("j", ""));
         double d = Double.parseDouble(parts[3].replace("k", ""));
         return new Quaternion(a, b, c, d);
      } catch (NumberFormatException e) {
         throw new IllegalArgumentException("Invalid input string: " + s);
      }
   }

   /** Clone of the quaternion.
    * @return independent clone of <code>this</code>
    */
   @Override
   public Object clone() throws CloneNotSupportedException {
      return new Quaternion(r, i, j, k);
   }

   /** Test whether the quaternion is zero.
    * @return true, if the real part and all the imaginary parts are (close to) zero
    */
   public boolean isZero() {
      return Math.abs(r) < EPS && Math.abs(i) < EPS && Math.abs(j) < EPS && Math.abs(k) < EPS;
   }

   /** Conjugate of the quaternion. Expressed by the formula
    *     conjugate(a+bi+cj+dk) = a-bi-cj-dk
    * @return conjugate of <code>this</code>
    */
   public Quaternion conjugate() {
      return new Quaternion(r, -i, -j, -k);
   }

   /** Opposite of the quaternion. Expressed by the formula
    *    opposite(a+bi+cj+dk) = -a-bi-cj-dk
    * @return quaternion <code>-this</code>
    */
   public Quaternion opposite() {
      return new Quaternion(-r, -i, -j, -k);
   }

   /** Sum of quaternions. Expressed by the formula
    *    (a1+b1i+c1j+d1k) + (a2+b2i+c2j+d2k) = (a1+a2) + (b1+b2)i + (c1+c2)j + (d1+d2)k
    * @param q addend
    * @return quaternion <code>this+q</code>
    */
   public Quaternion plus (Quaternion q) {
      return new Quaternion(r + q.r, i + q.i, j + q.j, k + q.k);
   }

   /** Product of quaternions. Expressed by the formula
    *  (a1+b1i+c1j+d1k) * (a2+b2i+c2j+d2k) = (a1a2-b1b2-c1c2-d1d2) + (a1b2+b1a2+c1d2-d1c2)i +
    *  (a1c2-b1d2+c1a2+d1b2)j + (a1d2+b1c2-c1b2+d1a2)k
    * @param q factor
    * @return quaternion <code>this*q</code>
    */
   public Quaternion times (Quaternion q) {
      double a = r * q.r - i * q.i - j * q.j - k * q.k;
      double b = r * q.i + i * q.r + j * q.k - k * q.j;
      double c = r * q.j - i * q.k + j * q.r + k * q.i;
      double d = r * q.k + i * q.j - j * q.i + k * q.r;
      return new Quaternion(a, b, c, d);
   }

   /** Multiplication by a coefficient.
    * @param r coefficient
    * @return quaternion <code>this*r</code>
    */
   public Quaternion times (double r) {
      double a = getRpart() * r;
      double b = getIpart() * r;
      double c = getJpart() * r;
      double d = getKpart() * r;

      return new Quaternion(a, b, c, d);
   }

   /** Inverse of the quaternion. Expressed by the formula
    *     1/(a+bi+cj+dk) = a/(a*a+b*b+c*c+d*d) +
    *     ((-b)/(a*a+b*b+c*c+d*d))i + ((-c)/(a*a+b*b+c*c+d*d))j + ((-d)/(a*a+b*b+c*c+d*d))k
    * @return quaternion <code>1/this</code>
    */
   public Quaternion inverse() {
      double normSquared = getRpart() * getRpart() + getIpart() * getIpart()
              + getJpart() * getJpart() + getKpart() * getKpart();
      if (isZero()) {
         throw new ArithmeticException("Quaternion is zero, cannot be inverted.");
      }
      return new Quaternion(getRpart() / normSquared, -getIpart() / normSquared,
              -getJpart() / normSquared, -getKpart() / normSquared);
   }

   /** Difference of quaternions. Expressed as addition to the opposite.
    * @param q subtrahend
    * @return quaternion <code>this-q</code>
    */
   public Quaternion minus (Quaternion q) {
      return new Quaternion(getRpart() - q.getRpart(), getIpart() - q.getIpart(),
              getJpart() - q.getJpart(), getKpart() - q.getKpart());
   }

   /** Right quotient of quaternions. Expressed as multiplication to the inverse.
    * @param q (right) divisor
    * @return quaternion <code>this*inverse(q)</code>
    */
   public Quaternion divideByRight (Quaternion q) {
      if (q.isZero()) {
         throw new ArithmeticException("Right operand is 0, cannot divide by zero.");
      }
      return times(q.inverse());
   }

   /** Left quotient of quaternions.
    * @param q (left) divisor
    * @return quaternion <code>inverse(q)*this</code>
    */
   public Quaternion divideByLeft (Quaternion q) {
      if (q.isZero()) {
         throw new ArithmeticException("Left operand is 0, cannot divide by zero.");
      }
      return q.inverse().times(this);
   }

   /** Equality test of quaternions. Difference of equal numbers
    *     is (close to) zero.
    * @param qo second quaternion
    * @return logical value of the expression <code>this.equals(qo)</code>
    */
   @Override
   public boolean equals (Object qo) {
      if (qo instanceof Quaternion) {
         Quaternion q = (Quaternion) qo;
         return Math.abs(getRpart() - q.getRpart()) < EPS &&
                 Math.abs(getIpart() - q.getIpart()) < EPS &&
                 Math.abs(getJpart() - q.getJpart()) < EPS &&
                 Math.abs(getKpart() - q.getKpart()) < EPS;
      }
      return false;
   }

   /** Dot product of quaternions. (p*conjugate(q) + q*conjugate(p))/2
    * @param q factor
    * @return dot product of this and q
    */
   public Quaternion dotMult (Quaternion q) {
      Quaternion pConj = conjugate();
      Quaternion qConj = q.conjugate();
      return pConj.times(q).plus(qConj.times(this)).times(0.5);
   }

   public Quaternion pow(int n){
      if (n == 0) {
         return new Quaternion(1, 0, 0, 0);
      } else if (n == 1) {
         return new Quaternion(r, i, j, k);
      } else if (n == -1) {
         return inverse();
      } else if (n > 1) {
         Quaternion q = new Quaternion(r, i, j, k);
         for (int i = 1; i < n; i++) {
            q = q.times(this);
         }
         return q;
      } else { // n < -1
         return pow(-n).inverse();
      }
   }

   /** Integer hashCode has to be the same for equal objects.
    * @return hashcode
    */
   @Override
   public int hashCode() {
      return Objects.hash(r, i, j, k);
   }

   /** Norm of the quaternion. Expressed by the formula
    *     norm(a+bi+cj+dk) = Math.sqrt(a*a+b*b+c*c+d*d)
    * @return norm of <code>this</code> (norm is a real number)
    */
   public double norm() {
      return Math.sqrt(r*r + i*i + j*j + k*k);
   }

   /** Main method for testing purposes.
    * @param args command line parameters
    */
   public static void main (String[] args) {
      Quaternion q1 = new Quaternion(1, 2, 3, 4);
      Quaternion q2 = new Quaternion(5, 6, 7, 8);

      System.out.println("q1 = " + q1);
      System.out.println("q2 = " + q2);

      Quaternion sum = q1.plus(q2);
      System.out.println("q1 + q2 = " + sum);

      Quaternion product = q1.times(q2);
      System.out.println("q1 * q2 = " + product);

      Quaternion inverse = q1.inverse();
      System.out.println("Inverse of q1 = " + inverse);

      Quaternion dotProduct = q1.dotMult(q2);
      System.out.println("Dot product of q1 and q2 = " + dotProduct);

      boolean equal = q1.equals(q2);
      System.out.println("q1 equals q2? " + equal);

      double norm = q1.norm();
      System.out.println("Norm of q1 = " + norm);
   }

}