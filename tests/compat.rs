//! Compatibility tests for `rsomics-plink-model` against PLINK 1.9 `--model`.
//!
//! The golden test runs our binary against a committed PLINK fileset and
//! diffs the output against PLINK's own `.model` (captured at fixture-build
//! time). The oracle test, when a `plink` binary is on PATH, runs both live
//! and requires byte-identical `.model` output.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-model"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn run_ours(prefix: &Path) -> String {
    let out = Command::new(bin())
        .arg(prefix)
        .output()
        .expect("run rsomics-plink-model");
    assert!(
        out.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8 output")
}

#[test]
fn golden_matches_committed_plink_output() {
    let prefix = golden_dir().join("sample");
    let ours = run_ours(&prefix);
    let expected = std::fs::read_to_string(golden_dir().join("sample.model.expected"))
        .expect("read expected .model");
    assert_eq!(
        ours, expected,
        "output diverged from committed PLINK .model"
    );
}

#[test]
fn matches_plink_oracle_when_available() {
    let plink = which_plink();
    let Some(plink) = plink else {
        eprintln!("SKIP: no `plink` binary on PATH or ~/oracle-bin — oracle compat skipped");
        return;
    };
    let ver = Command::new(&plink).arg("--version").output().unwrap();
    let ver = String::from_utf8_lossy(&ver.stdout);
    assert!(
        ver.contains("PLINK v1.9"),
        "oracle is not PLINK v1.9: {ver}"
    );

    let dir = tempfile::tempdir_in(scratch()).expect("tempdir");
    let prefix = golden_dir().join("sample");
    let out_prefix = dir.path().join("oracle");
    let status = Command::new(&plink)
        .args(["--bfile"])
        .arg(&prefix)
        .args(["--model", "--out"])
        .arg(&out_prefix)
        .status()
        .expect("run plink");
    assert!(status.success(), "plink --model failed");

    let oracle = std::fs::read_to_string(out_prefix.with_extension("model")).expect("read oracle");
    let ours = run_ours(&prefix);
    assert_eq!(ours, oracle, "rsomics-plink-model diverged from PLINK 1.9");
}

fn which_plink() -> Option<PathBuf> {
    let home = std::env::var("HOME").ok()?;
    let cand = PathBuf::from(&home).join("oracle-bin/plink");
    if cand.exists() {
        return Some(cand);
    }
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|p| p.join("plink"))
        .find(|p| p.exists())
}

fn scratch() -> PathBuf {
    std::env::var_os("TMPDIR")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
}
