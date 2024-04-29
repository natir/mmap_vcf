/* crate use */
use biotest;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use biotest::Format as _;

/* project use */
use biommap;

// Create reader
#[biommap::derive::sequential_parser(name = CountFastqRecord, data_type = u64, block_type = biommap::block::Block<memmap2::Mmap>, block_producer = biommap::fastq::File2Block, record_producer = biommap::fastq::Block2Record)]
fn parser(&mut self, _record: biommap::fastq::Record, counter: &mut u64) {
    *counter += 1;
}

pub fn block_size(c: &mut Criterion) {
    let mut g = c.benchmark_group("fastq");

    let mut rng = biotest::rand();
    let temp = tempfile::NamedTempFile::new().unwrap();
    let generator = biotest::Fastq::builder().sequence_len(150).build().unwrap();

    generator.create(temp.path(), &mut rng, 1000).unwrap();

    for pow in 10..18 {
        let block_size = 2u64.pow(pow);

        g.bench_with_input(
            BenchmarkId::new("block_size", block_size),
            &block_size,
            |b, block_size| {
                b.iter(|| {
                    let mut records_counter = 0;

                    let mut parser = CountFastqRecord::new();

                    parser
                        .with_blocksize(*block_size, temp.path(), &mut records_counter)
                        .unwrap();
                    black_box(records_counter);
                })
            },
        );
    }
}

criterion_group!(benches, block_size);
criterion_main!(benches);
