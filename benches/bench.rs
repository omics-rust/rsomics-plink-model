use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pgen::Pgen;
use rsomics_plink_model::{DEFAULT_CELL, model_test};
use std::hint::black_box;
use std::path::PathBuf;

fn bench_model(c: &mut Criterion) {
    let prefix = std::env::var("PLINK_MODEL_FIXTURE")
        .unwrap_or_else(|_| "/data3/liangjy/tmp/plink-model/perf".to_string());
    let pgen = Pgen::load(&PathBuf::from(&prefix)).expect("load fixture");
    c.bench_function("model", |b| {
        b.iter(|| black_box(model_test(black_box(&pgen), DEFAULT_CELL)));
    });
}

criterion_group!(benches, bench_model);
criterion_main!(benches);
