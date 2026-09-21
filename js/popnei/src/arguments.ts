/**
 * What the package checks of an argument before it reaches the core.
 *
 * A number of JavaScript is a float64, and the code wasm-bindgen generates
 * turns it into the integer the Rust takes by throwing its fraction away
 * and keeping it modulo 2^32, with no error: a ploidy of 2.5 and one of
 * 2^32 + 2 both arrive as 2, a size of 2^32 + 1 as blocks of one variant,
 * and -1 as 4294967295. Bytes that are not a `Uint8Array` are read as
 * whatever their memory holds. Python refuses all of them, so the package
 * refuses them here, before the call, and says what was given.
 */

/**
 * The largest number the package hands to the core, 2^32 - 1.
 *
 * A whole number of Rust is 32 bits wide in wasm, and what the generated
 * code does with a larger one is to keep it modulo 2^32: 2^32 + 2 would
 * arrive as a ploidy of 2, and 2^32 + 1 as blocks of one variant.
 */
const LARGEST_WHOLE_NUMBER = 4294967295;

/**
 * `value` when it is a whole number of 1 or more that the core holds, and
 * an `Error` that names `argument` and what was given otherwise.
 *
 * @throws {Error} When `value` is not such a number.
 */
export function wholeNumberOfOneOrMore(
  argument: string,
  value: unknown,
): number {
  if (
    typeof value !== "number" ||
    !Number.isSafeInteger(value) ||
    value < 1 ||
    value > LARGEST_WHOLE_NUMBER
  ) {
    throw new Error(
      `popnei: \`${argument}\` is a whole number of 1 or more and at most ` +
        `${LARGEST_WHOLE_NUMBER}, and ${whatWasGiven(value)} was given`,
    );
  }
  return value;
}

/**
 * `value` when it is a boolean, and an `Error` otherwise.
 *
 * @throws {Error} When `value` is not a boolean.
 */
export function aBoolean(argument: string, value: unknown): boolean {
  if (typeof value !== "boolean") {
    throw new Error(
      `popnei: \`${argument}\` is true or false, and ${whatWasGiven(value)} was given`,
    );
  }
  return value;
}

/**
 * The bytes of `value` when it is a `Uint8Array`, and an `Error` otherwise.
 *
 * @throws {Error} When `value` is not a `Uint8Array`.
 */
export function bytes(argument: string, value: unknown): Uint8Array {
  if (!(value instanceof Uint8Array)) {
    throw new Error(
      `popnei: \`${argument}\` is a Uint8Array with the bytes of the file, and ` +
        `${whatWasGiven(value)} was given; text is turned into bytes with ` +
        "new TextEncoder().encode(text), and a file of node is read with " +
        'new Uint8Array(await readFile(path))',
    );
  }
  return value;
}

/**
 * The names of `value` when it is an array of strings, and an `Error`
 * otherwise.
 *
 * One name written where the array goes, `fields: "alleles"`, is the case
 * this catches: a string spread into an array is its letters, and popnei
 * would look for a field called `a`.
 *
 * @throws {Error} When `value` is not an array of strings.
 */
export function namesOfFields(argument: string, value: unknown): string[] {
  if (!Array.isArray(value) || value.some((name) => typeof name !== "string")) {
    throw new Error(
      `popnei: \`${argument}\` is an array of names, and ${whatWasGiven(value)} ` +
        `was given; one field is asked for with ${argument}: ` +
        `${typeof value === "string" ? `["${value}"]` : '["chrom"]'}`,
    );
  }
  return value as string[];
}

/** What was given, for the message of an argument that was refused. */
export function whatWasGiven(value: unknown): string {
  if (typeof value === "string") {
    return `the string \`${value}\``;
  }
  if (value === null) {
    return "null";
  }
  if (typeof value === "object") {
    const name: unknown = (value as { constructor?: { name?: string } })
      .constructor?.name;
    return typeof name === "string" ? `a ${name}` : "an object";
  }
  return `the ${typeof value} ${String(value)}`;
}
