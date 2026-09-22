use kspike::*;
use std::time::Instant;
fn time<F: FnMut() -> (Vec<u32>, Vec<u32>)>(name: &str, reps: usize, mut f: F) -> (Vec<u32>, Vec<u32>) {
    let mut best = f64::MAX; let mut out = None;
    for _ in 0..reps { let t = Instant::now(); let r = f(); best = best.min(t.elapsed().as_secs_f64()); out = Some(r); }
    println!("  {name:<28} {best:.4} s"); out.unwrap()
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (nv, ni, na): (usize, usize, i8) = (a[1].parse().unwrap(), a[2].parse().unwrap(), a[3].parse().unwrap());
    let with_naive = a.get(4).map(|s| s == "naive").unwrap_or(false);
    println!("block of {nv} variants x {ni} individuals, {na} alleles, 3% missing");
    let g = make_block(nv, ni, na, 0.03, 42);
    let t = Instant::now(); let b = pack(&g, nv, ni); println!("  pack alone                   {:.4} s ({} sets)", t.elapsed().as_secs_f64(), b.nsets);
    let rb = time("bits, 1 thread", 3, || bits(&g, nv, ni));
    let rp = time("bits, rayon all cores", 3, || bits_par(&g, nv, ni));
    assert!(rb == rp);
    let rm = time("faer f32 products, 1 thread", 2, || matmul(&g, nv, ni));
    assert!(rb == rm, "matmul differs");
    if with_naive { let rn = time("pair by pair on i8, 1 thread", 1, || naive(&g, nv, ni)); assert!(rb == rn, "naive differs"); }
    println!("  checksum {}", checksum_pub(&rb));
}
