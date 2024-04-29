//! An efficient bioinformatics file parser based on memory mapping of file.

#![warn(missing_docs)]

/* std use */

/* crate use */

use clap::Parser as _;

/* project use */

/// Example: Count fasta record in file.
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "biommap",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    /// Input path
    #[clap(short = 'i', long = "input")]
    pub input_path: std::path::PathBuf,

    // Basic parameter
    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    pub verbosity: u8,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    pub ts: Option<stderrlog::Timestamp>,
}

fn main() -> anyhow::Result<()> {
    // parse cli
    let params = Command::parse();

    // Setup logger
    stderrlog::new()
        .quiet(params.quiet)
        .verbosity(params.verbosity as usize)
        .timestamp(params.ts.unwrap_or(stderrlog::Timestamp::Off))
        .init()
        .unwrap();

    let mut records_counter = 0;

    let mut reader = needletail::parse_fastx_file(&params.input_path)?;

    while let Some(result) = reader.next() {
        let _record = result?;

        records_counter += 1;
    }

    println!(
        "{} contains {} records",
        params.input_path.into_os_string().into_string().unwrap(),
        records_counter,
    );

    Ok(())
}
