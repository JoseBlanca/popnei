/**
 * `writeVars` and `writeVcf` with `onBytes`, which give the file to a
 * function of the application in pieces while the pass runs, as "The file
 * of a writer in pieces" of `docs/specs/js_sources.md` has it.
 *
 * What a call with the option gives is checked against the same call
 * without it: the pieces joined are its bytes and the counts of the pass
 * are its counts. The memory of wasm the option saves is measured in
 * `write_in_pieces_memory.test.ts`, alone in its file.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Progress, Variants, VarsWritten } from "popnei";
import { init, openVars, openVcf, writeVars, writeVcf } from "popnei";

import loadTheWasm from "../wasm/popnei.js";
import { referenceVcf, vcfOfDrawnGenotypes } from "./reference.ts";

await init();

/** The WebAssembly of the core, whose memory no piece may be a view into. */
const wasm = await loadTheWasm();

/** How many bytes every piece but the last holds, 1 MiB. */
const BYTES_PER_PIECE = 1048576;

/** The value the application throws to stop a writer. */
const THE_CANCEL = { cancelled: "by the test" };

/**
 * The VCF of 14000 variants of 600 individuals that `vars_memory.test.ts`
 * writes, 33.6 MB of text, which is several pieces written plain and as a
 * vars file.
 */
const DRAWN = vcfOfDrawnGenotypes(14000, 600);

/** The pieces a call gave, in their order, and those pieces joined. */
function piecesGivenTo(
  writes: (onBytes: (piece: Uint8Array) => void) => unknown,
): { pieces: Uint8Array[]; joined: Uint8Array; result: unknown } {
  const pieces: Uint8Array[] = [];
  const result = writes((piece) => {
    pieces.push(piece);
  });
  const joined = new Uint8Array(
    pieces.reduce((numBytes, piece) => numBytes + piece.length, 0),
  );
  let at = 0;
  for (const piece of pieces) {
    joined.set(piece, at);
    at += piece.length;
  }
  return { pieces, joined, result };
}

/** What `run` threw, or a failed assertion when it threw nothing. */
function whatWasThrownBy(run: () => unknown): unknown {
  try {
    run();
  } catch (thrown) {
    return thrown;
  }
  assert.fail("the call threw nothing");
}

/**
 * That every piece but the last holds 1 MiB and the last one from 1 byte to
 * 1 MiB.
 */
function assertTheSizesOf(pieces: readonly Uint8Array[]): void {
  assert.ok(pieces.length >= 1, "the function was never called");
  for (const piece of pieces.slice(0, -1)) {
    assert.equal(piece.length, BYTES_PER_PIECE);
  }
  const last = pieces.at(-1)?.length ?? 0;
  assert.ok(last >= 1 && last <= BYTES_PER_PIECE, `the last piece is ${last} bytes`);
}

/** The writers of the tests, each with the call it makes with and without the pieces. */
const THE_WRITES: readonly {
  name: string;
  whole: (variants: Variants) => { bytes: Uint8Array; passStats: unknown };
  inPieces: (
    variants: Variants,
    onBytes: (piece: Uint8Array) => void,
  ) => { passStats: unknown };
}[] = [
  {
    name: "writeVars",
    whole: (variants) => writeVars(variants, { numVarsPerBlock: 3 }),
    inPieces: (variants, onBytes) =>
      writeVars(variants, { numVarsPerBlock: 3, onBytes }),
  },
  {
    name: "writeVcf bgzipped",
    whole: (variants) => writeVcf(variants),
    inPieces: (variants, onBytes) => writeVcf(variants, { onBytes }),
  },
  {
    name: "writeVcf plain",
    whole: (variants) => writeVcf(variants, { bgzip: false }),
    inPieces: (variants, onBytes) =>
      writeVcf(variants, { bgzip: false, onBytes }),
  },
];

for (const write of THE_WRITES) {
  for (const file of ["cases.vcf", "many.vcf"]) {
    test(`${write.name} of ${file} gives in one piece the bytes it gives whole`, async () => {
      const variants = openVcf(await referenceVcf(file), { onlyPassed: false });
      try {
        const whole = write.whole(variants);
        const { pieces, joined, result } = piecesGivenTo((onBytes) =>
          write.inPieces(variants, onBytes),
        );
        assert.equal(pieces.length, 1);
        assert.deepEqual(joined, whole.bytes);
        assert.deepEqual(result, { passStats: whole.passStats });
      } finally {
        variants.free();
      }
    });
  }
}

for (const write of THE_WRITES) {
  test(`${write.name} of a vars file gives in one piece the bytes it gives whole`, async () => {
    const vcf = openVcf(await referenceVcf("many.vcf"), { onlyPassed: false });
    const variants = openVars(writeVars(vcf).bytes);
    vcf.free();
    try {
      const whole = write.whole(variants);
      const { pieces, joined, result } = piecesGivenTo((onBytes) =>
        write.inPieces(variants, onBytes),
      );
      assert.equal(pieces.length, 1);
      assert.deepEqual(joined, whole.bytes);
      assert.deepEqual(result, { passStats: whole.passStats });
    } finally {
      variants.free();
    }
  });
}

for (const write of THE_WRITES) {
  test(`${write.name} of 14000 variants gives the bytes it gives whole in pieces of 1 MiB`, () => {
    const variants = openVcf(DRAWN, { onlyPassed: false });
    try {
      const whole = write.whole(variants);
      const { pieces, joined, result } = piecesGivenTo((onBytes) =>
        write.inPieces(variants, onBytes),
      );
      assert.ok(pieces.length >= 3, `${pieces.length} pieces`);
      assertTheSizesOf(pieces);
      assert.deepEqual(joined, whole.bytes);
      assert.deepEqual(result, { passStats: whole.passStats });
    } finally {
      variants.free();
    }
  });
}

test("the pieces are the application's and stay as they were after a second write of the whole file", () => {
  const variants = openVcf(DRAWN, { onlyPassed: false });
  try {
    const whole = writeVcf(variants, { bgzip: false });
    const { pieces } = piecesGivenTo((onBytes) =>
      writeVcf(variants, { bgzip: false, onBytes }),
    );
    // The same pass again, whole, grows the memory of wasm by the file, and
    // a piece that were a view into it would have been moved or written over.
    writeVcf(variants, { bgzip: false });
    let at = 0;
    for (const piece of pieces) {
      assert.notEqual(piece.buffer, wasm.memory.buffer);
      assert.deepEqual(piece, whole.bytes.subarray(at, at + piece.length));
      at += piece.length;
    }
  } finally {
    variants.free();
  }
});

for (const write of THE_WRITES) {
  test(`${write.name} throws the value onBytes threw at the second piece, and the next call runs whole`, () => {
    const variants = openVcf(DRAWN, { onlyPassed: false });
    try {
      let calls = 0;
      const thrown = whatWasThrownBy(() =>
        write.inPieces(variants, () => {
          calls += 1;
          if (calls === 2) {
            throw THE_CANCEL;
          }
        }),
      );
      assert.equal(thrown, THE_CANCEL);
      assert.equal(calls, 2);
      const whole = write.whole(variants);
      const { joined } = piecesGivenTo((onBytes) =>
        write.inPieces(variants, onBytes),
      );
      assert.deepEqual(joined, whole.bytes);
    } finally {
      variants.free();
    }
  });
}

for (const write of THE_WRITES) {
  test(`free from inside onBytes while ${write.name} writes is refused`, async () => {
    const variants = openVcf(await referenceVcf("many.vcf"), { onlyPassed: false });
    try {
      const thrown = whatWasThrownBy(() =>
        write.inPieces(variants, () => {
          variants.free();
        }),
      );
      assert.ok(thrown instanceof Error, `it threw ${String(thrown)}`);
      assert.match(thrown.message, /a run is reading these variants/);
      assert.equal(write.whole(variants).bytes.length > 0, true);
    } finally {
      variants.free();
    }
  });
}

test("a VCF with a wrong line after 14000 good ones gives pieces and then the error of that line", () => {
  const wrong = new TextEncoder().encode("chr1\t14001\t.\tA\tC\t.\tPASS\t.\tGT\tnot_a_genotype\n");
  const vcf = new Uint8Array(DRAWN.length + wrong.length);
  vcf.set(DRAWN);
  vcf.set(wrong, DRAWN.length);
  const variants = openVcf(vcf, { onlyPassed: false });
  try {
    let calls = 0;
    const thrown = whatWasThrownBy(() =>
      writeVcf(variants, {
        bgzip: false,
        onBytes: () => {
          calls += 1;
        },
      }),
    );
    assert.ok(calls >= 1, "no piece was given before the wrong line");
    assert.ok(thrown instanceof Error, `it threw ${String(thrown)}`);
    // Two lines of header and 14000 variants come before it.
    assert.match(thrown.message, /^line 14003 of the VCF/);
  } finally {
    variants.free();
  }
});

for (const write of ["writeVars", "writeVcf"] as const) {
  test(`${write} refuses an onBytes that is not a function before the pass starts`, async () => {
    const variants = openVcf(await referenceVcf("many.vcf"));
    try {
      const told: Progress[] = [];
      variants.onProgress((progress) => {
        told.push(progress);
      });
      const writer = write === "writeVars" ? writeVars : writeVcf;
      assert.throws(
        // @ts-expect-error: the option is a function, and a page may give it
        // anything.
        () => writer(variants, { onBytes: 3 }),
        /`onBytes` is a function/,
      );
      assert.equal(told.length, 0);
    } finally {
      variants.free();
    }
  });
}

test("an onBytes of undefined is one not given", async () => {
  const variants = openVcf(await referenceVcf("many.vcf"));
  try {
    const written = writeVcf(variants, { onBytes: undefined });
    assert.ok(written.bytes.length > 0);
  } finally {
    variants.free();
  }
});

test("a file of exactly two pieces is given as two pieces of 1 MiB and no empty one", async () => {
  const many = new TextDecoder().decode(await referenceVcf("many.vcf"));
  const firstLineEnd = many.indexOf("\n") + 1;
  const padding = "##padding=";
  const numPadded =
    2 * BYTES_PER_PIECE - new TextEncoder().encode(many).length - padding.length - 1;
  assert.ok(numPadded > 0);
  const padded = new TextEncoder().encode(
    `${many.slice(0, firstLineEnd)}${padding}${"x".repeat(numPadded)}\n${many.slice(firstLineEnd)}`,
  );
  assert.equal(padded.length, 2 * BYTES_PER_PIECE);
  const variants = openVcf(padded, { onlyPassed: false });
  try {
    const { pieces, joined } = piecesGivenTo((onBytes) =>
      writeVcf(variants, { bgzip: false, onBytes }),
    );
    assert.deepEqual(
      pieces.map((piece) => piece.length),
      [BYTES_PER_PIECE, BYTES_PER_PIECE],
    );
    assert.deepEqual(joined, padded);
  } finally {
    variants.free();
  }
});

for (const write of THE_WRITES) {
  test(`${write.name} refuses an onBytes that returns a promise, after its first call`, async () => {
    const variants = openVcf(await referenceVcf("many.vcf"), { onlyPassed: false });
    try {
      let calls = 0;
      const thrown = whatWasThrownBy(() =>
        write.inPieces(variants, async () => {
          calls += 1;
        }),
      );
      assert.equal(calls, 1);
      assert.ok(thrown instanceof Error, `it threw ${String(thrown)}`);
      assert.match(thrown.message, /`onBytes`.*promise/);
      assert.ok(write.whole(variants).bytes.length > 0);
    } finally {
      variants.free();
    }
  });
}

test("a piece the heap has no room for throws the error of the browser and leaves the variants free to go", () => {
  const variants = openVcf(DRAWN, { onlyPassed: false });
  const TheUint8Array = globalThis.Uint8Array;
  const noRoom = new RangeError("Array buffer allocation failed");
  let numPieces = 0;
  globalThis.Uint8Array = new Proxy(TheUint8Array, {
    construct(target, args, newTarget) {
      // A piece is built from its length or from a view of its bytes.
      const first: unknown = args[0];
      const length =
        typeof first === "number" ? first : (first as { length?: number })?.length;
      if (length === BYTES_PER_PIECE) {
        numPieces += 1;
        if (numPieces >= 3) {
          throw noRoom;
        }
      }
      return Reflect.construct(target, args, newTarget);
    },
  });
  let thrown: unknown;
  try {
    thrown = whatWasThrownBy(() =>
      writeVcf(variants, { bgzip: false, onBytes: () => {} }),
    );
  } finally {
    globalThis.Uint8Array = TheUint8Array;
  }
  assert.equal(thrown, noRoom);
  variants.free();
});

test("writeVars and writeVcf run from inside onBytes over the same variants", () => {
  const variants = openVcf(DRAWN, { onlyPassed: false });
  try {
    const whole = writeVcf(variants, { bgzip: false });
    const vars = writeVars(variants).bytes;
    const inside: number[] = [];
    const { joined } = piecesGivenTo((onBytes) =>
      writeVcf(variants, {
        bgzip: false,
        onBytes: (piece) => {
          if (inside.length === 0) {
            inside.push(writeVars(variants).bytes.length);
            const { joined: innerVcf } = piecesGivenTo((innerOnBytes) =>
              writeVcf(variants, { bgzip: false, onBytes: innerOnBytes }),
            );
            inside.push(innerVcf.length);
          }
          onBytes(piece);
        },
      }),
    );
    assert.deepEqual(inside, [vars.length, whole.bytes.length]);
    assert.deepEqual(joined, whole.bytes);
  } finally {
    variants.free();
  }
});

test("the result of a call with onBytes has no bytes, and one without it is a VarsWritten", async () => {
  const variants = openVcf(await referenceVcf("many.vcf"), { onlyPassed: false });
  try {
    const inPieces = writeVars(variants, { onBytes: () => {} });
    // @ts-expect-error: a call with `onBytes` gives the counts alone.
    assert.equal(inPieces.bytes, undefined);
    const vcfInPieces = writeVcf(variants, { onBytes: () => {} });
    // @ts-expect-error: a call with `onBytes` gives the counts alone.
    assert.equal(vcfInPieces.bytes, undefined);
    const whole: VarsWritten = writeVars(variants);
    const notGiven: VarsWritten = writeVars(variants, { onBytes: undefined });
    assert.equal(whole.bytes.length, notGiven.bytes.length);
  } finally {
    variants.free();
  }
});
