// The draws of SplitMix64 that docs/specs/filters.md checks the filter that
// keeps variants at random against, from java.util.SplittableRandom, which is
// that generator. Run from the root of the repository, with OpenJDK 26.0.2.1,
// the version that gave the numbers of the spec:
//
//     java tests/reference/filters/SplitMix.java
//
// It prints the first five draws from a seed of 1234567, nextLong, and the
// first ten numbers from 0 to 1 from a seed of 42, nextDouble, which is the
// top 53 bits of a draw divided by 2^53, to six decimals.
import java.util.Locale;
import java.util.SplittableRandom;

public class SplitMix {
    public static void main(String[] args) {
        SplittableRandom from1234567 = new SplittableRandom(1234567L);
        for (int draw = 0; draw < 5; draw++) {
            System.out.println(Long.toUnsignedString(from1234567.nextLong()));
        }
        SplittableRandom from42 = new SplittableRandom(42L);
        StringBuilder numbers = new StringBuilder();
        for (int number = 0; number < 10; number++) {
            numbers.append(String.format(Locale.ROOT, "%.6f ", from42.nextDouble()));
        }
        System.out.println(numbers.toString().trim());
    }
}
