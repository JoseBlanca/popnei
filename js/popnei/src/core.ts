/**
 * What the two entry points of the package share: the functions of the
 * core crate, and the check that the WebAssembly is in place before one of
 * them is called.
 *
 * `node.ts` and `web.ts` differ only in where they get the bytes of the
 * WebAssembly from, and each of them calls `wasmIsLoaded` when it has
 * them. Neither the functions nor the check are exported to a user
 * directly: both entry points re-export what a user calls.
 */

import { version as versionOfTheCoreCrate } from "../wasm/popnei.js";

let theWasmIsLoaded = false;

/** Marks the WebAssembly of the core as loaded, called by `init`. */
export function wasmIsLoaded(): void {
  theWasmIsLoaded = true;
}

/**
 * The version of the core crate, `major.minor.patch`, which is the version
 * of this package.
 *
 * A user who reports a result names the code that gave it with this, so it
 * is the version of the Rust that did the work and not one written again
 * here.
 *
 * @throws {Error} When `init` has not been awaited.
 */
export function version(): string {
  if (!theWasmIsLoaded) {
    throw new Error(
      "popnei: await init() before calling any other function of the package",
    );
  }
  return versionOfTheCoreCrate();
}
