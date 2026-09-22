/**
 * What the package checks of an argument before it reaches the core.
 *
 * The bytes of a file are also asked about: the code wasm-bindgen
 * generates allocates the whole length of a `Uint8Array` inside the memory
 * of wasm before any code of popnei runs, and an allocation that fails
 * there is a trap, which ends the module. `room_for_bytes` of the binding
 * crate asks for that memory first and gives an `Error` when it is not
 * there.
 *
 * A number of JavaScript is a float64, and the code wasm-bindgen generates
 * turns it into the integer the Rust takes by throwing its fraction away
 * and keeping it modulo 2^32, with no error: a ploidy of 2.5 and one of
 * 2^32 + 2 both arrive as 2, a size of 2^32 + 1 as blocks of one variant,
 * and -1 as 4294967295. Bytes that are not a `Uint8Array` are read as
 * whatever their memory holds. Python refuses all of them, so the package
 * refuses them here, before the call, and says what was given.
 */

import { room_for_bytes as roomForBytes } from "../wasm/popnei.js";

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
 * `value` when it is a whole number of 0 or more that the core holds, and an
 * `Error` that names `argument` and what was given otherwise.
 *
 * What it refuses is what would reach the core as another number: the
 * generated code throws the fraction of 2.5 away, keeps 2^32 + 2 modulo
 * 2^32 and turns a negative number into a count of about four thousand
 * million. Whether the number is one the argument takes, a histogram of 1
 * bin or more and a ploidy of 1 to 255, is a rule of the core, which holds
 * for every pass and not for this call alone.
 *
 * The arguments that come through here all take 0: `numPrinComps`, how many
 * components the weights are asked for, where 0 asks for none;
 * `minNumSnps`, how many variants a pair of individuals needs before it
 * gets a distance, where 0 is every pair that was called at all and where
 * pyNei takes a negative number and does with it what it does with 0; and
 * `minNumIndividuals`, `ploidy` and `histKwargs.numBins` of the per variant
 * statistics, whose own rules the core holds.
 *
 * @throws {Error} When `value` is not such a number.
 */
export function wholeNumberOfZeroOrMore(
  argument: string,
  value: unknown,
): number {
  if (
    typeof value !== "number" ||
    !Number.isSafeInteger(value) ||
    value < 0 ||
    value > LARGEST_WHOLE_NUMBER
  ) {
    throw new Error(
      `popnei: \`${argument}\` is a whole number of 0 or more and at most ` +
        `${LARGEST_WHOLE_NUMBER}, and ${whatWasGiven(value)} was given`,
    );
  }
  return value;
}

/**
 * `value` when it is a string, and an `Error` that names `argument` and what
 * was given otherwise.
 *
 * Which strings the argument takes, the two kinds of bins of a histogram
 * among them, is a rule of the core: it refuses a name it does not know and
 * writes the ones it knows in the message. What the generated code does with
 * what is no string at all is to throw a `TypeError` of its own, which names
 * neither the argument nor what was given.
 *
 * @throws {Error} When `value` is not a string.
 */
export function aString(argument: string, value: unknown): string {
  if (typeof value !== "string") {
    throw new Error(
      `popnei: \`${argument}\` is a name, and ${whatWasGiven(value)} was given`,
    );
  }
  return value;
}

/**
 * `value` when it is a number, and an `Error` that names `argument` and what
 * was given otherwise.
 *
 * It is the threshold of a filter that comes through here, and what it
 * refuses is what is not a number at all: whether the number is one the
 * filter takes, from 0 to 1, is the rule of the core, which says it of the
 * threshold of every pass and not of this call alone. A call with no
 * threshold gives `undefined`, and the code wasm-bindgen generates would
 * hand the core a NaN for it, `null` as a threshold of 0 and the string
 * `"0.5"` as 0.5, each of them with no error.
 *
 * @throws {Error} When `value` is not a number.
 */
export function aNumber(argument: string, value: unknown): number {
  if (typeof value !== "number") {
    throw new Error(
      `popnei: \`${argument}\` is a number, and ${whatWasGiven(value)} was given`,
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
 * The bytes of `value` when it is a `Uint8Array` that can be read and that
 * the memory of wasm takes, and an `Error` otherwise.
 *
 * @throws {Error} When `value` is not a `Uint8Array`, when its buffer was
 * transferred, which leaves the array with nothing to read, and when the
 * memory of wasm does not take a copy of it.
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
  // A page that sends bytes to a web worker transfers their buffer, which
  // leaves the array it came from with a length of 0 and no memory behind
  // it. The code wasm-bindgen generates throws a `TypeError` of its own on
  // one, which names neither the argument nor what happened to it.
  // `detached` is of node 26 and of the browsers of 2024; where it is not
  // there, this passes and the array is read as the empty one it looks
  // like.
  if ((value.buffer as { detached?: unknown }).detached === true) {
    throw new Error(
      `popnei: the buffer of \`${argument}\` was transferred, to a web worker ` +
        "or somewhere else, and the bytes of the file are there and not in " +
        "this array; the worker that was given them is where they are read",
    );
  }
  // The copy into the memory of wasm is made by the generated code, before
  // any code of popnei runs, and an allocation that fails there is a trap
  // that leaves the module unusable. This asks for the memory first, and
  // what it grew is what that copy then finds.
  roomForBytes(value.length);
  return value;
}

/** What one of the names is and an example of one, for the message. */
export interface WhatTheNamesAre {
  /** What one name is, `field` or `individual`. */
  oneOfThem: string;
  /**
   * A name of that kind, `chrom` or `ind00`, which the message writes
   * inside the array the user should have written. What they gave is
   * written there instead when it is one name as a string.
   */
  anExample: string;
}

/**
 * The values of `value` when it is a `Float64Array` that holds a table of
 * `numRows` rows of `numCols` values each, row after row, and that the
 * memory of wasm takes, and an `Error` otherwise.
 *
 * The two sides are checked against the values there are here and not in the
 * core: a table given with fewer rows than it has would be the analysis of
 * the first of them, with nothing to show it, and one given with more traits
 * than it has would read the second row as the end of the first. The core
 * refuses only the second of the two, so a caller that says `numCols: 3` of
 * a table of 4 traits has to be stopped before the call.
 *
 * @throws {Error} When `value` is not a `Float64Array`, when its buffer was
 * transferred, which leaves the array with nothing to read, when it does not
 * hold `numRows` times `numCols` values, and when the memory of wasm does not
 * take a copy of it.
 */
export function valuesOfATable(
  argument: string,
  value: unknown,
  numRows: number,
  numCols: number,
): Float64Array {
  if (!(value instanceof Float64Array)) {
    throw new Error(
      `popnei: \`${argument}\` is a Float64Array with the values of the table, ` +
        `row after row, and ${whatWasGiven(value)} was given; an array of ` +
        `numbers is turned into one with Float64Array.from(numbers)`,
    );
  }
  // A page that sends the values to a web worker transfers their buffer,
  // which leaves the array it came from with a length of 0 and no memory
  // behind it, as `bytes` says above.
  if ((value.buffer as { detached?: unknown }).detached === true) {
    throw new Error(
      `popnei: the buffer of \`${argument}\` was transferred, to a web worker ` +
        "or somewhere else, and the values of the table are there and not in " +
        "this array; the worker that was given them is where they are read",
    );
  }
  if (value.length !== numRows * numCols) {
    throw new Error(
      `popnei: the table was given as ${numRows} x ${numCols}, which is ` +
        `${numRows * numCols} values, and \`${argument}\` holds ${value.length}`,
    );
  }
  // The copy into the memory of wasm is made by the generated code, before
  // any code of popnei runs, and an allocation that fails there is a trap
  // that leaves the module unusable, as `bytes` says above. A value is 8
  // bytes.
  roomForBytes(value.length * 8);
  return value;
}

/**
 * The names of `value` when it is an array of strings, and an `Error`
 * otherwise, which says what one of them is, a `field` or an `individual`,
 * and writes an example of it in what the user should have written.
 *
 * One name written where the array goes, `fields: "alleles"`, is the case
 * this catches: a string spread into an array is its letters, and popnei
 * would look for a field called `a`.
 *
 * @throws {Error} When `value` is not an array of strings.
 */
export function namesOf(
  argument: string,
  value: unknown,
  whatTheNamesAre: WhatTheNamesAre,
): string[] {
  if (!Array.isArray(value) || value.some((name) => typeof name !== "string")) {
    const { oneOfThem, anExample } = whatTheNamesAre;
    throw new Error(
      `popnei: \`${argument}\` is an array of names, and ${whatWasGiven(value)} ` +
        `was given; one ${oneOfThem} is asked for with ${argument}: ` +
        `["${typeof value === "string" ? value : anExample}"]`,
    );
  }
  return value as string[];
}

/**
 * What was given, for the message of an argument that was refused.
 *
 * The name of the class of an object goes inside the words `of the type`,
 * which is where the Python package puts it too: the article that would
 * come before it is `a` for a `Uint8Array` and `an` for an `Array`, and no
 * rule of the letters tells the two apart.
 */
export function whatWasGiven(value: unknown): string {
  if (typeof value === "string") {
    return `the string \`${value}\``;
  }
  if (value === null) {
    return "null";
  }
  if (value === undefined) {
    return "undefined";
  }
  if (typeof value === "object") {
    const name: unknown = (value as { constructor?: { name?: string } })
      .constructor?.name;
    return typeof name === "string"
      ? `an object of the type \`${name}\``
      : "an object of no type";
  }
  return `the ${typeof value} ${String(value)}`;
}
