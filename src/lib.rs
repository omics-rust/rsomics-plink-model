//! PLINK1 `--model` full case/control genotypic association test.
//!
//! Emits five per-variant tests (GENO, TREND, ALLELIC, DOM, REC) matching
//! `plink --model`.

#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]

pub mod model;
pub mod stats;

pub use model::{DEFAULT_CELL, VariantTests, model_test, print_model};
pub use stats::chi2_sf;
