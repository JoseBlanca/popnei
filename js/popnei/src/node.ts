/**
 * The entry point under node, where the WebAssembly of the core is read
 * from the file beside this package.
 *
 * node resolves the package to this file through the condition `node` of
 * the field `exports` of `package.json`; a bundler and a page get `web.ts`
 * instead. The two export the same functions.
 */

import { readFile } from "node:fs/promises";

import loadTheWasm from "../wasm/popnei.js";
import { wasmIsLoaded } from "./core.js";

export type { Block, Field } from "./block.js";
export { version } from "./core.js";
export { calcPairwiseKosmanDists, Distances } from "./dists.js";
export type { CalcPairwiseKosmanDistsOptions } from "./dists.js";
export type { FilteringStats, Step } from "./filters.js";
export { calcGwas } from "./gwas.js";
export type {
  CalcGwasOptions,
  GwasModel,
  GwasNullModel,
  GwasResult,
  GwasStats,
  TestType,
  TraitType,
} from "./gwas.js";
export { calcKinship, Kinship } from "./kinship.js";
export type { CalcKinshipOptions, KinshipPcsResult } from "./kinship.js";
export { calcLdAndDistPerPop, calcRogersHuffR2Matrix } from "./ld.js";
export type {
  CalcLdAndDistPerPopOptions,
  CalcRogersHuffR2MatrixOptions,
  LdAndDistPerPop,
  LdBins,
  R2Matrix,
} from "./ld.js";
export { openVars, writeVars } from "./io_vars.js";
export type { VarsWritten, WriteVarsOptions } from "./io_vars.js";
export { openVcf } from "./io_vcf.js";
export type { OpenVcfOptions } from "./io_vcf.js";
export { doPca, doPcaFromVariants } from "./pca.js";
export type {
  DoPcaFromVariantsOptions,
  DoPcaOptions,
  PcaResult,
  VariantsPcaResult,
} from "./pca.js";
export { calcPopDists } from "./pop_dists.js";
export type {
  CalcPopDistsOptions,
  PopDistGroup,
  PopDistMeasure,
  PopDists,
} from "./pop_dists.js";
export { calcPerIndividualStats, calcPerVarDistribs } from "./stats.js";
export type {
  BinType,
  HistKwargs,
  PerIndividualStats,
  PerVarDistribs,
  PerVarDistribsOptions,
  PerVarStat,
  PolyVarsStats,
  StatsDistrib,
} from "./stats.js";
export { Variants } from "./variant.js";
export type { Blocks, IterBlocksOptions, PassStats } from "./variant.js";

let loading: Promise<void> | undefined;

/**
 * Loads the WebAssembly of the core, and returns when it can be called.
 *
 * It has to be awaited before any other function of the package, and a
 * second call returns the same promise as the first, so that the
 * WebAssembly is loaded once.
 */
export function init(): Promise<void> {
  loading ??= readFile(new URL("../wasm/popnei_bg.wasm", import.meta.url))
    .then((wasm) => loadTheWasm({ module_or_path: wasm }))
    .then(wasmIsLoaded);
  return loading;
}
