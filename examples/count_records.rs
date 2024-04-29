//! An efficient bioinformatics file parser based on memory mapping of file.

#![warn(missing_docs)]

/* std use */

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use clap::Parser as _;

/* project use */

/// Enum to select sequence type
#[derive(Default, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum SequenceType {
    /// Fasta sequence type
    #[default]
    Fasta,

    /// Fastq sequence type
    Fastq,
}

impl std::fmt::Display for SequenceType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SequenceType::Fasta => write!(f, "fasta"),
            SequenceType::Fastq => write!(f, "fastq"),
        }
    }
}

#[cfg(not(feature = "parallel"))]
#[biommap::derive::sequential_parser(name = CountFastaRecord, data_type = u64, block_type = biommap::block::Block<memmap2::Mmap>, block_producer = biommap::fasta::File2Block, record_producer = biommap::fasta::Block2Record)]
fn parser(&mut self, _record: biommap::fasta::Record, counter: &mut u64) {
    *counter += 1;
}

#[cfg(feature = "parallel")]
#[biommap::derive::sharedstate_parser(name = CountFastaRecord, data_type = std::sync::atomic::AtomicU64, block_producer = biommap::fasta::File2Block, record_producer = biommap::fasta::Block2Record)]
fn parser(_record: biommap::fasta::Record, counter: &std::sync::atomic::AtomicU64) {
    counter.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
}

#[cfg(not(feature = "parallel"))]
#[biommap::derive::sequential_parser(name = CountFastqRecord, data_type = u64, block_type = biommap::block::Block<memmap2::Mmap>, block_producer = biommap::fastq::File2Block, record_producer = biommap::fastq::Block2Record)]
fn parser(&mut self, _record: biommap::fastq::Record, counter: &mut u64) {
    *counter += 1;
}

#[cfg(feature = "parallel")]
#[biommap::derive::sharedstate_parser(name = CountFastqRecord, data_type = std::sync::atomic::AtomicU64, block_producer = biommap::fastq::File2Block, record_producer = biommap::fastq::Block2Record)]
fn parser(_record: biommap::fastq::Record, counter: &std::sync::atomic::AtomicU64) {
    counter.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
}

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

    /// Block size
    #[clap(short = 'b', long = "block-size")]
    pub block_size: u64,

    /// Sequence type
    #[clap(short = 's', long = "type", value_enum)]
    pub sequence_type: SequenceType,

    #[cfg(feature = "parallel")]
    /// Number of thread usable
    #[clap(short = 't', long = "threads")]
    pub threads: usize,

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

    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads)
        .build_global()?;

    #[cfg(feature = "parallel")]
    let records_counter = std::sync::atomic::AtomicU64::new(0);

    #[cfg(not(feature = "parallel"))]
    let mut records_counter = 0;

    match params.sequence_type {
        SequenceType::Fasta => {
            log::info!("fasta mode");

            let mut parser = CountFastaRecord::new();

            #[cfg(feature = "parallel")]
            parser.with_blocksize(params.block_size, &params.input_path, &records_counter)?;

            #[cfg(not(feature = "parallel"))]
            parser.with_blocksize(params.block_size, &params.input_path, &mut records_counter)?;
        }
        SequenceType::Fastq => {
            log::info!("fastq mode");

            let mut parser = CountFastqRecord::new();

            #[cfg(feature = "parallel")]
            parser.with_blocksize(params.block_size, &params.input_path, &records_counter)?;

            #[cfg(not(feature = "parallel"))]
            parser.with_blocksize(params.block_size, &params.input_path, &mut records_counter)?;
        }
    }

    #[cfg(feature = "parallel")]
    let value = records_counter.load(std::sync::atomic::Ordering::SeqCst);
    #[cfg(not(feature = "parallel"))]
    let value = records_counter;

    println!(
        "{} contains {} {} records",
        params.input_path.into_os_string().into_string().unwrap(),
        value,
        params.sequence_type,
    );

    Ok(())
}
