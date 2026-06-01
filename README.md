# rsomics-plink-model

Full case/control **genotypic** association test — a Rust port of PLINK 1.9
`--model`. For every variant it builds the 2×3 case/control × genotype table
and emits five tests:

| TEST | Table | DF | Statistic |
|---|---|---|---|
| `GENO` | 2×3 genotype (HomA1 / Het / HomA2) | 2 | Pearson χ² |
| `TREND` | allele counts, genotype scores 0/1/2 | 1 | Cochran-Armitage trend |
| `ALLELIC` | 2×2 allele | 1 | Pearson χ² |
| `DOM` | (HomA1+Het) vs HomA2 | 1 | Pearson χ² |
| `REC` | HomA1 vs (Het+HomA2) | 1 | Pearson χ² |

`GENO`, `DOM` and `REC` report `NA` when any cell of the 2×3 genotype table
is below the `--cell` threshold (default 5, matching PLINK); `TREND` and
`ALLELIC` are always computed. This is the distinct genotypic-model test, not
the single allelic χ² of `rsomics-plink-assoc` (`--assoc`).

## Usage

```sh
rsomics-plink-model PREFIX            # reads PREFIX.{bed,bim,fam}
rsomics-plink-model PREFIX --cell 5   # cell-count threshold (default 5)
```

Output is the PLINK `.model` table on stdout:

```
 CHR  SNP   A1   A2     TEST            AFF          UNAFF        CHISQ   DF            P
   1  rs3    A    G     GENO       26/87/54       15/87/65        3.968    2       0.1375
   1  rs3    A    G    TREND        139/195        117/217        3.414    1      0.06466
   1  rs3    A    G  ALLELIC        139/195        117/217        3.065    1      0.07998
   1  rs3    A    G      DOM         113/54         102/65         1.58    1       0.2088
   1  rs3    A    G      REC         26/141         15/152        3.364    1      0.06663
```

Phenotype comes from the `.fam` sixth column: `2` = affected (case), `1` =
unaffected (control); any other value excludes the sample. `AFF`/`UNAFF` hold
the genotype counts (`GENO`) or collapsed allele/genotype counts per test.

## Origin

This crate is an independent Rust reimplementation of PLINK 1.9 `--model`
based on:

- The published method: Purcell et al. 2007 (PLINK,
  [doi:10.1086/519795](https://doi.org/10.1086/519795)) and Chang et al. 2015
  (PLINK 1.9, [doi:10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8)).
- The Cochran-Armitage trend test (Armitage 1955).
- The public PLINK 1.9 binary file-format spec
  (<https://www.cog-genomics.org/plink/1.9/formats>) and the `--model`
  documentation (<https://www.cog-genomics.org/plink/1.9/assoc>).
- Black-box behaviour testing against the PLINK 1.9 binary
  (v1.9.0-b.7.7), including exact `.model` column layout, χ² formulae, the
  Cochran-Armitage statistic, the `--cell` NA threshold, and PLINK's `%g`-style
  number formatting.

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated from a fixed seed.

License: MIT OR Apache-2.0.
Upstream credit: PLINK 1.9 (Christopher Chang et al., GPLv3),
<https://www.cog-genomics.org/plink/1.9/>.
