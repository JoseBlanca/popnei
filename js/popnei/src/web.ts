/**
 * The entry point in a browser and in a bundler, where the WebAssembly of
 * the core is fetched from the address it sits at beside the JavaScript
 * that wasm-bindgen generated.
 *
 * A bundler resolves the package to this file through the field `exports`
 * of `package.json`, node gets `node.ts` instead, and the two export the
 * same functions. Whether the wasm file ends up where that loader looks
 * for it is the bundler's doing, and they differ: the README of the
 * package has what vite and esbuild were seen to do with the
 * `new URL("popnei_bg.wasm", import.meta.url)` that the loader holds.
 */

import loadTheWasm from "../wasm/popnei.js";
import { wasmIsLoaded } from "./core.js";

export type { Block, Field } from "./block.js";
export { version } from "./core.js";
export type { FilteringStats, Step } from "./filters.js";
export { openVars, writeVars } from "./io_vars.js";
export type { VarsWritten, WriteVarsOptions } from "./io_vars.js";
export { openVcf } from "./io_vcf.js";
export type { OpenVcfOptions } from "./io_vcf.js";
export { Variants } from "./variant.js";
export type { Blocks, IterBlocksOptions, PassStats } from "./variant.js";

let loading: Promise<void> | undefined;

/**
 * Fetches the WebAssembly of the core, and returns when it can be called.
 *
 * It has to be awaited before any other function of the package, and a
 * second call returns the same promise as the first, so that the
 * WebAssembly is fetched once.
 */
export function init(): Promise<void> {
  loading ??= loadTheWasm().then(wasmIsLoaded);
  return loading;
}
