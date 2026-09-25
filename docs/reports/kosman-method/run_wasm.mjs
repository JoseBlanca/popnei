import fs from 'node:fs';
const bytes = fs.readFileSync(process.argv[2]);
const { instance } = await WebAssembly.instantiate(bytes, {});
const { run, run_make } = instance.exports;
for (const [nv, ni, na] of [[5000, 1000, 2], [5000, 1000, 4], [500, 5000, 2]]) {
  let t = performance.now(); run_make(nv, ni, na); const make = (performance.now() - t) / 1000;
  for (const [m, name] of [[1, 'bits'], [2, 'faer f32 products']]) {
    let best = 1e9, cs;
    for (let r = 0; r < 2; r++) { t = performance.now(); cs = run(m, nv, ni, na); best = Math.min(best, (performance.now() - t) / 1000 - make); }
    console.log(`wasm, ${nv} x ${ni}, ${na} alleles, ${name}: ${best.toFixed(4)} s  (checksum ${cs})`);
  }
}
