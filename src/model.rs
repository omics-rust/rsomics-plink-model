//! PLINK1 `--model` full genotypic association.
//!
//! Per variant, build the 2×3 case/control × {HomA1, Het, HomA2} genotype
//! table, then emit five tests: GENO (2df genotypic Pearson), TREND
//! (Cochran-Armitage, weights 0/1/2), ALLELIC (2×2 allele Pearson), DOM
//! ({HomA1+Het} vs HomA2), REC (HomA1 vs {Het+HomA2}). DF=2 for GENO, else 1.
//!
//! Phenotype: .fam col 6, "2" = case (AFF), "1" = control (UNAFF); other =
//! excluded. A1 is the first .bim allele. GENO/DOM/REC report NA when any of
//! the six genotype-table cells is below the cell-count threshold (PLINK
//! `--cell`, default 5); TREND and ALLELIC are always computed.

use crate::stats::chi2_sf;
use rsomics_pgen::{Pgen, PgenMmap, Sample, Variant};
use std::io::Write;

pub const DEFAULT_CELL: u32 = 5;

/// Variant-major genotype access shared by the in-memory and mmap readers, so
/// `--model` runs identically whichever one main feeds it.
pub trait BedRows {
    fn samples(&self) -> &[Sample];
    fn variants(&self) -> &[Variant];
    fn n_variants(&self) -> usize;
    fn variant_row(&self, v: usize) -> &[u8];
}

impl BedRows for Pgen {
    fn samples(&self) -> &[Sample] {
        &self.samples
    }
    fn variants(&self) -> &[Variant] {
        &self.variants
    }
    fn n_variants(&self) -> usize {
        self.n_variants()
    }
    fn variant_row(&self, v: usize) -> &[u8] {
        self.variant_row(v)
    }
}

impl BedRows for PgenMmap {
    fn samples(&self) -> &[Sample] {
        &self.samples
    }
    fn variants(&self) -> &[Variant] {
        &self.variants
    }
    fn n_variants(&self) -> usize {
        self.n_variants()
    }
    fn variant_row(&self, v: usize) -> &[u8] {
        self.variant_row(v)
    }
}

/// Counts of HomA1 / Het / HomA2 for one phenotype group.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
struct GenoCounts {
    hom_a1: u32,
    het: u32,
    hom_a2: u32,
}

/// The five `--model` tests in PLINK output order.
const TESTS: [&str; 5] = ["GENO", "TREND", "ALLELIC", "DOM", "REC"];

/// All per-variant results, borrowing the variant identifiers from the source
/// reader. The genotype counts are kept oriented (A1 = minor) so the count
/// columns can be formatted at print time; `a1`/`a2` are owned only when the
/// minor-allele swap relabels them. Only the five (chisq, df, p) statistics are
/// precomputed in the parallel build (`None` = NA), keeping it allocation-free.
pub struct VariantTests<'a> {
    pub chrom: &'a str,
    pub snp: &'a str,
    pub a1: &'a str,
    pub a2: &'a str,
    aff: GenoCounts,
    unaff: GenoCounts,
    stats: [Option<(f64, u32, f64)>; 5],
}

/// Per-sample phenotype class derived once from .fam.
#[derive(Clone, Copy)]
enum Pheno {
    Aff,
    Unaff,
    Skip,
}

/// Low bit set for every 2-bit lane: marks the lane positions whose low bit
/// carries a class/group indicator after the bit-plane split.
const LO_LANES: u64 = 0x5555_5555_5555_5555;

/// Phenotype group as a bit-plane mask aligned to the .bed 2-bit lane layout:
/// for each group, the low bit of a sample's lane is set when the sample is in
/// that group. Samples per 64-bit word: 32.
struct GroupMasks {
    case: Vec<u64>,
    ctrl: Vec<u64>,
}

impl GroupMasks {
    fn build(pheno: &[Pheno]) -> Self {
        let n_words = pheno.len().div_ceil(32);
        let mut case = vec![0u64; n_words];
        let mut ctrl = vec![0u64; n_words];
        for (s, &p) in pheno.iter().enumerate() {
            let bit = 1u64 << ((s % 32) * 2);
            match p {
                Pheno::Aff => case[s / 32] |= bit,
                Pheno::Unaff => ctrl[s / 32] |= bit,
                Pheno::Skip => {}
            }
        }
        Self { case, ctrl }
    }
}

/// Six masked class×group bit-streams per word, in `[case, ctrl] × [HomA1, Het,
/// HomA2]` order. Each carries one set bit per matching sample at that sample's
/// lane-low-bit position.
#[inline]
fn class_streams(word: u64, cm: u64, um: u64) -> [u64; 6] {
    let lo = word & LO_LANES;
    let hi = (word >> 1) & LO_LANES;
    let hom_a1 = !hi & !lo & LO_LANES;
    let het = hi & !lo;
    let hom_a2 = hi & lo;
    [
        hom_a1 & cm,
        het & cm,
        hom_a2 & cm,
        hom_a1 & um,
        het & um,
        hom_a2 & um,
    ]
}

/// Count, per phenotype group, the HomA1/Het/HomA2 genotypes in one variant
/// row, one 32-sample word at a time. The six class×group bit-streams reduce to
/// a `popcnt` each; missing lanes never enter a stream, so they are ignored.
fn count_row(row: &[u8], masks: &GroupMasks) -> (GenoCounts, GenoCounts) {
    let mut t = [0u32; 6];
    let mut w = 0;
    let chunks = row.chunks_exact(8);
    let tail = chunks.remainder();
    for chunk in chunks {
        let s = class_streams(
            u64::from_le_bytes(chunk.try_into().unwrap()),
            masks.case[w],
            masks.ctrl[w],
        );
        for k in 0..6 {
            t[k] += s[k].count_ones();
        }
        w += 1;
    }
    if !tail.is_empty() {
        let mut buf = [0u8; 8];
        buf[..tail.len()].copy_from_slice(tail);
        let s = class_streams(u64::from_le_bytes(buf), masks.case[w], masks.ctrl[w]);
        for k in 0..6 {
            t[k] += s[k].count_ones();
        }
    }
    (
        GenoCounts {
            hom_a1: t[0],
            het: t[1],
            hom_a2: t[2],
        },
        GenoCounts {
            hom_a1: t[3],
            het: t[4],
            hom_a2: t[5],
        },
    )
}

fn pheno_masks(pgen: &impl BedRows) -> GroupMasks {
    let pheno: Vec<Pheno> = pgen
        .samples()
        .iter()
        .map(|s| match s.phen.as_str() {
            "2" => Pheno::Aff,
            "1" => Pheno::Unaff,
            _ => Pheno::Skip,
        })
        .collect();
    GroupMasks::build(&pheno)
}

pub fn model_test<R: BedRows + Sync>(pgen: &R, cell: u32) -> Vec<VariantTests<'_>> {
    use rayon::prelude::*;
    let masks = pheno_masks(pgen);
    (0..pgen.n_variants())
        .into_par_iter()
        .map(|v| {
            let (aff, unaff) = count_row(pgen.variant_row(v), &masks);
            build_tests(&pgen.variants()[v], aff, unaff, cell)
        })
        .collect()
}

/// Compute and emit the full `.model` report in one fused streaming pass: each
/// block of variants is counted, tested, and rendered to text in parallel, then
/// written, so the 150 k×5-row table never materialises as an intermediate
/// `Vec<VariantTests>` nor gets a second cold-memory formatting pass.
pub fn write_model<R: BedRows + Sync>(
    pgen: &R,
    cell: u32,
    out: &mut impl Write,
) -> std::io::Result<()> {
    use rayon::prelude::*;
    let masks = pheno_masks(pgen);
    let snp_w = pgen
        .variants()
        .iter()
        .map(|v| v.id.len())
        .max()
        .map_or(4, |m| if m <= 4 { 4 } else { m + 2 });

    out.write_all(&header_line(snp_w))?;

    let line_slot = snp_w + 96;
    let block_slot = line_slot * 5;
    const BLOCK_VARIANTS: usize = 8192;
    let n = pgen.n_variants();
    let mut buf = vec![0u8; block_slot * n.min(BLOCK_VARIANTS)];

    for start in (0..n).step_by(BLOCK_VARIANTS) {
        let width = (start + BLOCK_VARIANTS).min(n) - start;
        let ends: Vec<usize> = buf
            .par_chunks_mut(block_slot)
            .take(width)
            .enumerate()
            .map(|(i, slot)| {
                let v = start + i;
                let (aff, unaff) = count_row(pgen.variant_row(v), &masks);
                let t = build_tests(&pgen.variants()[v], aff, unaff, cell);
                format_variant(slot, &t, snp_w)
            })
            .collect();
        let total: usize = ends.iter().sum();
        pack_blocks(&mut buf, block_slot, &ends, total);
        out.write_all(&buf[..total])?;
    }
    Ok(())
}

fn header_line(snp_w: usize) -> Vec<u8> {
    use std::io::Write as _;
    let mut header = Vec::with_capacity(96);
    writeln!(
        header,
        " CHR {:>snp_w$} {:>4} {:>4} {:>8} {:>14} {:>14} {:>12} {:>4} {:>12}",
        "SNP", "A1", "A2", "TEST", "AFF", "UNAFF", "CHISQ", "DF", "P"
    )
    .unwrap();
    header
}

impl GenoCounts {
    fn swap_homs(self) -> Self {
        Self {
            hom_a1: self.hom_a2,
            het: self.het,
            hom_a2: self.hom_a1,
        }
    }
}

fn build_tests(
    var: &rsomics_pgen::Variant,
    mut aff: GenoCounts,
    mut unaff: GenoCounts,
    cell: u32,
) -> VariantTests<'_> {
    // PLINK reports A1 as the minor allele. The .bed homozygous codes are tied
    // to the .bim allele order, so when the .bim A1 is the major allele we swap
    // the labels and the HomA1/HomA2 counts; an exact tie keeps .bim order.
    let bim_a1 = aff.hom_a1 + unaff.hom_a1;
    let bim_a2 = aff.hom_a2 + unaff.hom_a2;
    let (a1, a2) = if bim_a1 > bim_a2 {
        aff = aff.swap_homs();
        unaff = unaff.swap_homs();
        (var.a2.as_str(), var.a1.as_str())
    } else {
        (var.a1.as_str(), var.a2.as_str())
    };

    let geno_cells = [
        aff.hom_a1,
        aff.het,
        aff.hom_a2,
        unaff.hom_a1,
        unaff.het,
        unaff.hom_a2,
    ];
    let geno_ok = geno_cells.iter().all(|&c| c >= cell);

    let (aa1, aa2) = (aff.hom_a1 * 2 + aff.het, aff.het + aff.hom_a2 * 2);
    let (ua1, ua2) = (unaff.hom_a1 * 2 + unaff.het, unaff.het + unaff.hom_a2 * 2);

    let geno = geno_ok.then(|| {
        chi2_2x3(
            [aff.hom_a1, aff.het, aff.hom_a2],
            [unaff.hom_a1, unaff.het, unaff.hom_a2],
        )
    });
    let dom = geno_ok.then(|| {
        chi2_2x2(
            aff.hom_a1 + aff.het,
            aff.hom_a2,
            unaff.hom_a1 + unaff.het,
            unaff.hom_a2,
        )
    });
    let rec = geno_ok.then(|| {
        chi2_2x2(
            aff.hom_a1,
            aff.het + aff.hom_a2,
            unaff.hom_a1,
            unaff.het + unaff.hom_a2,
        )
    });

    let stat = |chisq: f64, df: u32| chisq.is_finite().then(|| (chisq, df, chi2_sf(chisq, df)));
    let stats = [
        geno.and_then(|c| stat(c, 2)),
        stat(cochran_armitage(aff, unaff), 1),
        stat(chi2_2x2(aa1, aa2, ua1, ua2), 1),
        dom.and_then(|c| stat(c, 1)),
        rec.and_then(|c| stat(c, 1)),
    ];

    VariantTests {
        chrom: var.chrom.as_str(),
        snp: var.id.as_str(),
        a1,
        a2,
        aff,
        unaff,
        stats,
    }
}

/// Pearson chi-squared on a 2×3 genotype table (df=2, no continuity correction).
fn chi2_2x3(aff: [u32; 3], unaff: [u32; 3]) -> f64 {
    let col: [f64; 3] = std::array::from_fn(|i| f64::from(aff[i] + unaff[i]));
    let n_aff: f64 = aff.iter().map(|&c| f64::from(c)).sum();
    let n_unaff: f64 = unaff.iter().map(|&c| f64::from(c)).sum();
    let n = n_aff + n_unaff;
    if n == 0.0 {
        return f64::NAN;
    }
    let mut chi = 0.0;
    for i in 0..3 {
        let e_aff = col[i] * n_aff / n;
        let e_unaff = col[i] * n_unaff / n;
        if e_aff > 0.0 {
            chi += (f64::from(aff[i]) - e_aff).powi(2) / e_aff;
        }
        if e_unaff > 0.0 {
            chi += (f64::from(unaff[i]) - e_unaff).powi(2) / e_unaff;
        }
    }
    chi
}

/// Pearson chi-squared on a 2×2 table (df=1, no continuity correction).
fn chi2_2x2(a: u32, b: u32, c: u32, d: u32) -> f64 {
    let (a, b, c, d) = (f64::from(a), f64::from(b), f64::from(c), f64::from(d));
    let n = a + b + c + d;
    let (r1, r2, c1, c2) = (a + b, c + d, a + c, b + d);
    if r1 == 0.0 || r2 == 0.0 || c1 == 0.0 || c2 == 0.0 {
        return f64::NAN;
    }
    let det = a * d - b * c;
    n * det * det / (r1 * r2 * c1 * c2)
}

/// Cochran-Armitage trend test with genotype scores (0, 1, 2), df=1.
/// T = N·[N·Σ wᵢrᵢ − R·Σ wᵢnᵢ]² / (R·S·[N·Σ wᵢ²nᵢ − (Σ wᵢnᵢ)²]).
fn cochran_armitage(aff: GenoCounts, unaff: GenoCounts) -> f64 {
    let r = [aff.hom_a1, aff.het, aff.hom_a2];
    let s = [unaff.hom_a1, unaff.het, unaff.hom_a2];
    let n_i: [f64; 3] = std::array::from_fn(|i| f64::from(r[i] + s[i]));
    let rr: f64 = r.iter().map(|&c| f64::from(c)).sum();
    let ss: f64 = s.iter().map(|&c| f64::from(c)).sum();
    let nn = rr + ss;
    let w = [0.0, 1.0, 2.0];
    let wr: f64 = (0..3).map(|i| w[i] * f64::from(r[i])).sum();
    let wn: f64 = (0..3).map(|i| w[i] * n_i[i]).sum();
    let w2n: f64 = (0..3).map(|i| w[i] * w[i] * n_i[i]).sum();
    let num = nn * wr - rr * wn;
    let den = rr * ss * (nn * w2n - wn * wn);
    if den == 0.0 {
        return f64::NAN;
    }
    nn * num * num / den
}

/// PLINK `.model` output; NA for non-finite/below-cell tests. The SNP column
/// is right-justified to a width that fits the longest variant ID (the same
/// rule PLINK uses): 4 for IDs up to 4 chars, otherwise the longest ID + 2.
///
/// Each variant's five-line block is formatted in parallel into its own
/// fixed-width slot; a gather pass packs the live bytes down to one contiguous
/// image that is emitted in a single write so the streamed bytes match PLINK
/// exactly while keeping the write path off the per-variant serial loop.
pub fn print_model(records: &[VariantTests], out: &mut impl Write) -> std::io::Result<()> {
    use rayon::prelude::*;
    let snp_w = snp_field_width(records);
    out.write_all(&header_line(snp_w))?;

    // The five fixed-width fields plus separators and newline never exceed this
    // for one line; the slack covers alleles or counts that overflow their pad.
    let line_slot = snp_w + 96;
    let block_slot = line_slot * 5;

    // Variants are blocked so a wide matrix never buffers the whole image at
    // once; each block's slots are filled in parallel, packed, then written.
    const BLOCK_VARIANTS: usize = 8192;
    let mut buf = vec![0u8; block_slot * records.len().min(BLOCK_VARIANTS)];

    for chunk in records.chunks(BLOCK_VARIANTS) {
        let ends: Vec<usize> = buf
            .par_chunks_mut(block_slot)
            .zip(chunk.par_iter())
            .map(|(slot, r)| format_variant(slot, r, snp_w))
            .collect();
        let total: usize = ends.iter().sum();
        pack_blocks(&mut buf, block_slot, &ends, total);
        out.write_all(&buf[..total])?;
    }
    Ok(())
}

/// Compact the per-variant slots (each `block_slot` wide, the i-th holding
/// `ends[i]` live bytes) down into the contiguous prefix `buf[..total]`.
fn pack_blocks(buf: &mut [u8], block_slot: usize, ends: &[usize], total: usize) {
    let mut write = 0usize;
    let mut read = 0usize;
    for &end in ends {
        if write != read {
            buf.copy_within(read..read + end, write);
        }
        write += end;
        read += block_slot;
    }
    debug_assert_eq!(write, total);
}

/// Format one variant's five-line block into the start of `slot`, returning the
/// number of bytes written.
fn format_variant(slot: &mut [u8], r: &VariantTests, snp_w: usize) -> usize {
    let (a, u) = (r.aff, r.unaff);
    let aa = [a.hom_a1 * 2 + a.het, a.het + a.hom_a2 * 2];
    let ua = [u.hom_a1 * 2 + u.het, u.het + u.hom_a2 * 2];
    let mut w = Slot { buf: slot, pos: 0 };
    for (i, &test) in TESTS.iter().enumerate() {
        w.padded(r.chrom.as_bytes(), 4);
        w.space();
        w.padded(r.snp.as_bytes(), snp_w);
        w.space();
        w.padded(r.a1.as_bytes(), 4);
        w.space();
        w.padded(r.a2.as_bytes(), 4);
        w.space();
        w.padded(test.as_bytes(), 8);
        w.space();

        let (aff_n, unaff_n): (&[u32], &[u32]) = match i {
            0 => (&[a.hom_a1, a.het, a.hom_a2], &[u.hom_a1, u.het, u.hom_a2]),
            1 | 2 => (&aa, &ua),
            3 => (&[a.hom_a1 + a.het, a.hom_a2], &[u.hom_a1 + u.het, u.hom_a2]),
            _ => (&[a.hom_a1, a.het + a.hom_a2], &[u.hom_a1, u.het + u.hom_a2]),
        };
        w.padded_slashed(aff_n, 14);
        w.space();
        w.padded_slashed(unaff_n, 14);
        w.space();

        match r.stats[i] {
            Some((chisq, df, p)) => {
                w.padded_sig(chisq, 12);
                w.space();
                w.padded_uint(u64::from(df), 4);
                w.space();
                w.padded_sig(p, 12);
            }
            None => {
                w.padded(b"NA", 12);
                w.space();
                w.padded(b"NA", 4);
                w.space();
                w.padded(b"NA", 12);
            }
        }
        w.newline();
    }
    w.pos
}

/// A cursor writing right-justified, space-separated fields straight into a
/// variant's pre-sized output slot, no intermediate per-field heap buffers.
struct Slot<'a> {
    buf: &'a mut [u8],
    pos: usize,
}

impl Slot<'_> {
    #[inline]
    fn space(&mut self) {
        self.buf[self.pos] = b' ';
        self.pos += 1;
    }

    #[inline]
    fn newline(&mut self) {
        self.buf[self.pos] = b'\n';
        self.pos += 1;
    }

    #[inline]
    fn pad(&mut self, content_len: usize, width: usize) {
        for _ in 0..width.saturating_sub(content_len) {
            self.buf[self.pos] = b' ';
            self.pos += 1;
        }
    }

    #[inline]
    fn raw(&mut self, s: &[u8]) {
        self.buf[self.pos..self.pos + s.len()].copy_from_slice(s);
        self.pos += s.len();
    }

    #[inline]
    fn padded(&mut self, s: &[u8], width: usize) {
        self.pad(s.len(), width);
        self.raw(s);
    }

    fn padded_uint(&mut self, n: u64, width: usize) {
        let mut tmp = [0u8; 20];
        let s = fmt_uint(&mut tmp, n);
        self.padded(s, width);
    }

    fn padded_slashed(&mut self, ns: &[u32], width: usize) {
        // Slash-joined counts never exceed a 2×3-table cell's worst case;
        // 32 bytes covers three 10-digit numbers plus two separators.
        let mut tmp = [0u8; 32];
        let mut len = 0;
        let mut digits = [0u8; 20];
        for (i, &n) in ns.iter().enumerate() {
            if i > 0 {
                tmp[len] = b'/';
                len += 1;
            }
            let d = fmt_uint(&mut digits, u64::from(n));
            tmp[len..len + d.len()].copy_from_slice(d);
            len += d.len();
        }
        self.padded(&tmp[..len], width);
    }

    fn padded_sig(&mut self, x: f64, width: usize) {
        let mut tmp = [0u8; 24];
        let len = fmt_sig(&mut tmp, x, 4);
        self.padded(&tmp[..len], width);
    }
}

fn snp_field_width(records: &[VariantTests]) -> usize {
    let longest = records.iter().map(|r| r.snp.len()).max().unwrap_or(0);
    if longest <= 4 { 4 } else { longest + 2 }
}

/// Significant digits of `x` (round-half-to-even on the shortest round-trip
/// decimal) and the base-10 exponent of the leading digit, written into a
/// fixed stack buffer to avoid per-value heap allocation.
struct SigDigits {
    digits: [u8; 17],
    len: usize,
    exp: i32,
}

/// Parse a ryu shortest-decimal string (`12.005`, `7.15e-5`, `1e-15`, `100.0`)
/// into significant digits (no leading or trailing zeros) and the base-10
/// exponent of the leading digit. `x` is non-zero and finite.
fn parse_shortest(s: &[u8]) -> ([u8; 17], usize, i32) {
    let mut raw = [0u8; 24];
    let mut raw_len = 0;
    let mut point = None; // index in `raw` just past the integer part
    let mut e_exp = 0i32;
    let mut i = 0;
    while i < s.len() {
        match s[i] {
            b'.' => point = Some(raw_len),
            b'e' => {
                let mut sign = 1i32;
                i += 1;
                if s[i] == b'-' {
                    sign = -1;
                    i += 1;
                } else if s[i] == b'+' {
                    i += 1;
                }
                let mut v = 0i32;
                while i < s.len() {
                    v = v * 10 + i32::from(s[i] - b'0');
                    i += 1;
                }
                e_exp = sign * v;
                break;
            }
            d => {
                raw[raw_len] = d;
                raw_len += 1;
            }
        }
        i += 1;
    }
    let int_digits = point.unwrap_or(raw_len);
    // Exponent of the first raw digit before trimming leading zeros.
    let mut exp = int_digits as i32 - 1 + e_exp;

    let mut start = 0;
    while start < raw_len - 1 && raw[start] == b'0' {
        start += 1;
        exp -= 1;
    }
    let mut end = raw_len;
    while end > start + 1 && raw[end - 1] == b'0' {
        end -= 1;
    }

    let mut digits = [0u8; 17];
    let len = end - start;
    digits[..len].copy_from_slice(&raw[start..end]);
    (digits, len, exp)
}

fn round_significant(x: f64, sig: usize) -> SigDigits {
    // ryu emits the shortest round-trip decimal (the same property PLINK's
    // dtoa relies on) far faster than the std Grisu path; we parse its output
    // into significant digits + a leading-digit exponent, then re-round.
    let mut ryu = ryu::Buffer::new();
    let s = ryu.format_finite(x.abs());
    let (mut digits, mut len, mut exp) = parse_shortest(s.as_bytes());

    if len > sig {
        let round_up = match digits[sig].cmp(&b'5') {
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => {
                digits[sig + 1..len].iter().any(|&d| d != b'0') || (digits[sig - 1] - b'0') % 2 == 1
            }
        };
        len = sig;
        if round_up {
            let mut i = sig;
            loop {
                if i == 0 {
                    digits.copy_within(0..sig, 1);
                    digits[0] = b'1';
                    exp += 1;
                    break;
                }
                i -= 1;
                if digits[i] == b'9' {
                    digits[i] = b'0';
                } else {
                    digits[i] += 1;
                    break;
                }
            }
        }
    }
    while len < sig {
        digits[len] = b'0';
        len += 1;
    }
    SigDigits { digits, len, exp }
}

/// Format `x` to `sig` significant digits PLINK-style (fixed-point for decimal
/// exponents in -4..sig, scientific otherwise; trailing zeros trimmed) into the
/// start of `out`, returning the number of bytes written. `out` must hold a
/// fully formatted value (24 bytes is ample for 4 significant digits).
fn fmt_sig(out: &mut [u8], x: f64, sig: usize) -> usize {
    let mut w = Cur { buf: out, pos: 0 };
    if !x.is_finite() {
        w.put_slice(b"NA");
        return w.pos;
    }
    if x == 0.0 {
        w.put(b'0');
        return w.pos;
    }
    let SigDigits { digits, len, exp } = round_significant(x, sig);

    if exp < -4 || exp >= sig as i32 {
        w.put(digits[0]);
        let frac_end = (1..len).rfind(|&i| digits[i] != b'0').map_or(0, |i| i + 1);
        if frac_end > 1 {
            w.put(b'.');
            w.put_slice(&digits[1..frac_end]);
        }
        w.put(b'e');
        w.put(if exp < 0 { b'-' } else { b'+' });
        let e = exp.unsigned_abs();
        if e < 10 {
            w.put(b'0');
        }
        let mut tmp = [0u8; 20];
        w.put_slice(fmt_uint(&mut tmp, u64::from(e)));
        return w.pos;
    }

    let point = exp + 1;
    if point <= 0 {
        w.put_slice(b"0.");
        w.put_zeros((-point) as usize);
        let frac_end = (0..len).rfind(|&i| digits[i] != b'0').map_or(0, |i| i + 1);
        w.put_slice(&digits[..frac_end]);
    } else if point as usize >= len {
        w.put_slice(&digits[..len]);
        w.put_zeros(point as usize - len);
    } else {
        let p = point as usize;
        w.put_slice(&digits[..p]);
        let frac_end = (p..len).rfind(|&i| digits[i] != b'0').map_or(p, |i| i + 1);
        if frac_end > p {
            w.put(b'.');
            w.put_slice(&digits[p..frac_end]);
        }
    }
    w.pos
}

/// A bare append cursor over a byte slice for the decimal formatter.
struct Cur<'a> {
    buf: &'a mut [u8],
    pos: usize,
}

impl Cur<'_> {
    #[inline]
    fn put(&mut self, b: u8) {
        self.buf[self.pos] = b;
        self.pos += 1;
    }
    #[inline]
    fn put_slice(&mut self, s: &[u8]) {
        self.buf[self.pos..self.pos + s.len()].copy_from_slice(s);
        self.pos += s.len();
    }
    #[inline]
    fn put_zeros(&mut self, n: usize) {
        self.buf[self.pos..self.pos + n].fill(b'0');
        self.pos += n;
    }
}

/// Decimal-format `n` into the back of `buf`, returning the live digit slice.
fn fmt_uint(buf: &mut [u8; 20], mut n: u64) -> &[u8] {
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    &buf[i..]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 5e-3, "got {a}, want {b}");
    }

    fn fmt_g(x: f64) -> String {
        let mut s = [0u8; 24];
        let len = fmt_sig(&mut s, x, 4);
        String::from_utf8(s[..len].to_vec()).unwrap()
    }

    /// Reference genotype count: one lane at a time, no carry-save tree.
    fn count_naive(row: &[u8], masks: &GroupMasks) -> (GenoCounts, GenoCounts) {
        let mut t = [0u32; 6];
        for (w, chunk) in row.chunks(8).enumerate() {
            let mut buf = [0u8; 8];
            buf[..chunk.len()].copy_from_slice(chunk);
            let s = class_streams(u64::from_le_bytes(buf), masks.case[w], masks.ctrl[w]);
            for k in 0..6 {
                t[k] += s[k].count_ones();
            }
        }
        let g = |i| GenoCounts {
            hom_a1: t[i],
            het: t[i + 1],
            hom_a2: t[i + 2],
        };
        (g(0), g(3))
    }

    #[test]
    fn count_row_matches_naive_across_word_counts() {
        // The golden fixture's rows are a handful of words; sweep odd sizes that
        // hit the zero-padded tail and many-word cases against an independent
        // chunk-at-a-time reference.
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut rng = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for n_samples in [1usize, 31, 32, 33, 256, 257, 1000, 8000, 8003] {
            let pheno: Vec<Pheno> = (0..n_samples)
                .map(|i| match i % 3 {
                    0 => Pheno::Aff,
                    1 => Pheno::Unaff,
                    _ => Pheno::Skip,
                })
                .collect();
            let masks = GroupMasks::build(&pheno);
            let row: Vec<u8> = (0..n_samples.div_ceil(4)).map(|_| rng() as u8).collect();
            assert_eq!(
                count_row(&row, &masks),
                count_naive(&row, &masks),
                "mismatch at n_samples={n_samples}"
            );
        }
    }

    #[test]
    fn statistics_match_plink() {
        // PLINK 1.9 --model output for a variant with case 26/87/54, control 15/87/65.
        let aff = GenoCounts {
            hom_a1: 26,
            het: 87,
            hom_a2: 54,
        };
        let unaff = GenoCounts {
            hom_a1: 15,
            het: 87,
            hom_a2: 65,
        };
        close(chi2_2x3([26, 87, 54], [15, 87, 65]), 3.968);
        close(cochran_armitage(aff, unaff), 3.414);
        close(chi2_2x2(139, 195, 117, 217), 3.065); // ALLELIC
        close(chi2_2x2(113, 54, 102, 65), 1.58); // DOM
        close(chi2_2x2(26, 141, 15, 152), 3.364); // REC
    }

    #[test]
    fn cell_threshold_gates_geno_dom_rec() {
        // A zero cell forces GENO/DOM/REC to NA; TREND/ALLELIC still emitted.
        let var = rsomics_pgen::Variant {
            chrom: "1".into(),
            id: "rs1".into(),
            cm: 0.0,
            pos: 1,
            a1: "A".into(),
            a2: "G".into(),
        };
        let aff = GenoCounts {
            hom_a1: 0,
            het: 15,
            hom_a2: 185,
        };
        let unaff = GenoCounts {
            hom_a1: 0,
            het: 13,
            hom_a2: 187,
        };
        let t = build_tests(&var, aff, unaff, DEFAULT_CELL);
        let na: Vec<&str> = TESTS
            .iter()
            .zip(t.stats)
            .filter(|(_, s)| s.is_none())
            .map(|(&name, _)| name)
            .collect();
        assert_eq!(na, ["GENO", "DOM", "REC"]);
    }

    #[test]
    fn snp_width_matches_plink() {
        let mk = |id: &'static str| VariantTests {
            chrom: "1",
            snp: id,
            a1: "A",
            a2: "G",
            aff: GenoCounts::default(),
            unaff: GenoCounts::default(),
            stats: [None; 5],
        };
        assert_eq!(snp_field_width(&[mk("rs1")]), 4);
        assert_eq!(snp_field_width(&[mk("abcd")]), 4);
        assert_eq!(snp_field_width(&[mk("abcde")]), 7);
        assert_eq!(snp_field_width(&[mk("rs100000")]), 10);
    }

    #[test]
    fn fmt_matches_plink_g_style() {
        assert_eq!(fmt_g(3.968), "3.968");
        assert_eq!(fmt_g(4.0), "4");
        assert_eq!(fmt_g(1.58), "1.58");
        assert_eq!(fmt_g(0.2918), "0.2918");
        assert_eq!(fmt_g(0.0006405), "0.0006405");
        assert_eq!(fmt_g(7.15e-05), "7.15e-05");
        assert_eq!(fmt_g(3.697e-08), "3.697e-08");
        assert_eq!(fmt_g(1.0), "1");
        // Round-half-to-even at exact .xxx5 decimal ties (PLINK dtoa_g parity).
        assert_eq!(fmt_g(12.005), "12");
        assert_eq!(fmt_g(22.445), "22.44");
        assert_eq!(fmt_g(65.025), "65.02");
        assert_eq!(fmt_g(85.805), "85.8");
        // Single-significant-digit scientific mantissa drops the decimal point.
        assert_eq!(fmt_g(8.0e-05), "8e-05");
        assert_eq!(fmt_g(1.0e-15), "1e-15");
    }
}
