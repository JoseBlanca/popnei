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

import type { Progress, Variants } from "popnei";
import { init, openVars, openVcf, writeVars, writeVcf } from "popnei";

import { referenceVcf, vcfOfDrawnGenotypes } from "./reference.ts";

await init();

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

test("writeVcf of a vars file gives in one piece the bytes it gives whole", async () => {
  const vcf = openVcf(await referenceVcf("many.vcf"), { onlyPassed: false });
  const variants = openVars(writeVars(vcf).bytes);
  vcf.free();
  try {
    const whole = writeVcf(variants, { bgzip: false });
    const { pieces, joined, result } = piecesGivenTo((onBytes) =>
      writeVcf(variants, { bgzip: false, onBytes }),
    );
    assert.equal(pieces.length, 1);
    assert.deepEqual(joined, whole.bytes);
    assert.deepEqual(result, { passStats: whole.passStats });
  } finally {
    variants.free();
  }
});

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

test("the pieces are the application's and stay as they were while the memory of wasm grows", () => {
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
      assert.ok(!(piece.buffer instanceof SharedArrayBuffer));
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
