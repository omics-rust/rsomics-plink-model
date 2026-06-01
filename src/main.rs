use clap::Parser;
use rsomics_pgen::Pgen;
use rsomics_plink_model::{DEFAULT_CELL, model_test, print_model};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "rsomics-plink-model",
    about = "PLINK1 full case/control genotypic association test (--model): GENO, TREND, ALLELIC, DOM, REC",
    version
)]
struct Cli {
    /// Path prefix for .bed/.bim/.fam (without extension)
    bfile: PathBuf,

    /// Minimum count per genotype-table cell below which GENO/DOM/REC report NA
    #[arg(long, default_value_t = DEFAULT_CELL)]
    cell: u32,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    let pgen = Pgen::load(&cli.bfile)?;
    let records = model_test(&pgen, cli.cell);
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    print_model(&records, &mut out)?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_definition_is_valid() {
        <Cli as clap::CommandFactory>::command().debug_assert();
    }
}
