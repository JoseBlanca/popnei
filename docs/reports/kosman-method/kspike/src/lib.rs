//! Spike: Kosman distances of every pair of individuals over one block,
//! three ways. gts is variants x individuals x 2, i8, -1 missing.

pub struct Rng(u64);
impl Rng {
    pub fn new(s: u64) -> Self { Rng(s) }
    pub fn next(&mut self) -> u64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; self.0 }
    pub fn unif(&mut self) -> f64 { (self.next() >> 11) as f64 / (1u64 << 53) as f64 }
}

/// Random block: each variant has an allele frequency, `num_alleles` alleles, `missing` rate per genotype.
pub fn make_block(nv: usize, ni: usize, num_alleles: i8, missing: f64, seed: u64) -> Vec<i8> {
    let mut r = Rng::new(seed);
    let mut g = vec![0i8; nv * ni * 2];
    for v in 0..nv {
        let p = 0.05 + 0.9 * r.unif();
        for i in 0..ni {
            let o = (v * ni + i) * 2;
            if r.unif() < missing { g[o] = -1; g[o + 1] = -1; continue; }
            for k in 0..2 {
                g[o + k] = if r.unif() < p { 0 } else { 1 + (r.next() % (num_alleles as u64 - 1)) as i8 };
            }
        }
    }
    g
}

fn npairs(ni: usize) -> usize { ni * (ni - 1) / 2 }

/// A: pair by pair on the int8, variants in the inner loop. Returns (2*dist_sum, n) per pair.
pub fn naive(g: &[i8], nv: usize, ni: usize) -> (Vec<u32>, Vec<u32>) {
    // transpose to individuals x variants x 2 so the inner loop is contiguous
    let mut t = vec![0i8; g.len()];
    for v in 0..nv { for i in 0..ni { let s = (v * ni + i) * 2; let d = (i * nv + v) * 2; t[d] = g[s]; t[d + 1] = g[s + 1]; } }
    let mut twice = vec![0u32; npairs(ni)]; let mut n = vec![0u32; npairs(ni)];
    let mut idx = 0;
    for i in 0..ni { let a = &t[i * nv * 2..(i + 1) * nv * 2];
        for j in i + 1..ni { let b = &t[j * nv * 2..(j + 1) * nv * 2];
            let (mut tw, mut nn) = (0u32, 0u32);
            for v in 0..nv {
                let (a0, a1, b0, b1) = (a[2 * v], a[2 * v + 1], b[2 * v], b[2 * v + 1]);
                let called = (a0 >= 0) & (a1 >= 0) & (b0 >= 0) & (b1 >= 0);
                let same = ((a0 == b0) & (a1 == b1)) | ((a0 == b1) & (a1 == b0));
                let share = (a0 == b0) | (a0 == b1) | (a1 == b0) | (a1 == b1);
                let d = if same { 0 } else if share { 1 } else { 2 };
                tw += if called { d } else { 0 }; nn += called as u32;
            }
            twice[idx] = tw; n[idx] = nn; idx += 1;
        }
    }
    (twice, n)
}

/// The bit sets of one block: per individual, `words` u64 per set; sets are called, then carries_a and hom_a per allele.
pub struct Bits { pub words: usize, pub nsets: usize, pub data: Vec<u64> }

pub fn pack(g: &[i8], nv: usize, ni: usize) -> Bits {
    let max_allele = g.iter().copied().max().unwrap_or(-1);
    let na = (max_allele as i32 + 1).max(0) as usize;
    let words = nv.div_ceil(64); let nsets = 1 + 2 * na;
    let mut data = vec![0u64; ni * nsets * words];
    for v in 0..nv { let (w, bit) = (v / 64, 1u64 << (v % 64));
        for i in 0..ni { let o = (v * ni + i) * 2; let (a0, a1) = (g[o], g[o + 1]);
            if a0 < 0 || a1 < 0 { continue; }
            let base = i * nsets * words;
            data[base + w] |= bit;
            data[base + (1 + 2 * a0 as usize) * words + w] |= bit;
            data[base + (1 + 2 * a1 as usize) * words + w] |= bit;
            if a0 == a1 { data[base + (2 + 2 * a0 as usize) * words + w] |= bit; }
        }
    }
    Bits { words, nsets, data }
}

#[inline]
fn pair(b: &Bits, i: usize, j: usize) -> (u32, u32) {
    let w = b.words; let stride = b.nsets * w;
    let a = &b.data[i * stride..(i + 1) * stride]; let c = &b.data[j * stride..(j + 1) * stride];
    let mut n = 0u32; for k in 0..w { n += (a[k] & c[k]).count_ones(); }
    let mut s = 0u32; for k in w..stride { s += (a[k] & c[k]).count_ones(); }
    (2 * n - s, n)
}

/// B: bit sets and popcount, one thread.
pub fn bits(g: &[i8], nv: usize, ni: usize) -> (Vec<u32>, Vec<u32>) {
    let b = pack(g, nv, ni);
    let mut twice = Vec::with_capacity(npairs(ni)); let mut n = Vec::with_capacity(npairs(ni));
    for i in 0..ni { for j in i + 1..ni { let (t, nn) = pair(&b, i, j); twice.push(t); n.push(nn); } }
    (twice, n)
}

#[cfg(not(target_family = "wasm"))]
pub fn bits_par(g: &[i8], nv: usize, ni: usize) -> (Vec<u32>, Vec<u32>) {
    use rayon::prelude::*;
    let b = pack(g, nv, ni);
    let rows: Vec<Vec<(u32, u32)>> = (0..ni).into_par_iter().map(|i| (i + 1..ni).map(|j| pair(&b, i, j)).collect()).collect();
    let mut twice = Vec::with_capacity(npairs(ni)); let mut n = Vec::with_capacity(npairs(ni));
    for r in rows { for (t, nn) in r { twice.push(t); n.push(nn); } }
    (twice, n)
}

/// C: pyNei's way, indicator matrices in f32 and products, with faer.
pub fn matmul(g: &[i8], nv: usize, ni: usize) -> (Vec<u32>, Vec<u32>) {
    use faer::Mat;
    let max_allele = g.iter().copied().max().unwrap_or(-1);
    let called = Mat::<f32>::from_fn(nv, ni, |v, i| { let o = (v * ni + i) * 2; (g[o] >= 0 && g[o + 1] >= 0) as u8 as f32 });
    let nmat = called.transpose() * &called;
    let mut acc = Mat::<f32>::zeros(ni, ni);
    for a in 0..=max_allele {
        let car = Mat::<f32>::from_fn(nv, ni, |v, i| { let o = (v * ni + i) * 2; (g[o] >= 0 && g[o + 1] >= 0 && (g[o] == a || g[o + 1] == a)) as u8 as f32 });
        acc += car.transpose() * &car;
        let hom = Mat::<f32>::from_fn(nv, ni, |v, i| { let o = (v * ni + i) * 2; (g[o] == a && g[o + 1] == a) as u8 as f32 });
        acc += hom.transpose() * &hom;
    }
    let mut twice = Vec::with_capacity(npairs(ni)); let mut n = Vec::with_capacity(npairs(ni));
    for i in 0..ni { for j in i + 1..ni { n.push(nmat[(i, j)] as u32); twice.push((2.0 * nmat[(i, j)] - acc[(i, j)]) as u32); } }
    (twice, n)
}

fn checksum(r: &(Vec<u32>, Vec<u32>)) -> f64 { r.0.iter().zip(&r.1).map(|(t, n)| *t as f64 / (2.0 * (*n).max(1) as f64)).sum() }

/// For wasm: method 0 naive, 1 bits, 2 matmul. Returns the checksum.
#[unsafe(no_mangle)]
pub extern "C" fn run(method: u32, nv: u32, ni: u32, num_alleles: u32) -> f64 {
    let (nv, ni) = (nv as usize, ni as usize);
    let g = make_block(nv, ni, num_alleles as i8, 0.03, 42);
    let r = match method { 0 => naive(&g, nv, ni), 1 => bits(&g, nv, ni), _ => matmul(&g, nv, ni) };
    checksum(&r)
}
/// Only the making of the block, to subtract.
#[unsafe(no_mangle)]
pub extern "C" fn run_make(nv: u32, ni: u32, num_alleles: u32) -> f64 { make_block(nv as usize, ni as usize, num_alleles as i8, 0.03, 42).len() as f64 }

pub fn checksum_pub(r: &(Vec<u32>, Vec<u32>)) -> f64 { checksum(r) }
