<h1 style="text-align: center;">biommap</h1>

![Test](https://github.com/natir/biommap/workflows/Test/badge.svg)
![Lints](https://github.com/natir/biommap/workflows/Lints/badge.svg)
![MSRV](https://github.com/natir/biommap/workflows/MSRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/biommap/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/biommap)
[![Documentation](https://github.com/natir/biommap/workflows/Documentation/badge.svg)](https://natir.github.io/biommap/biommap)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/biommap/blob/master/LICENSE)


An efficient bioinformatics file parser based on memory mapping of file.

**WARNING**:
- biommap work only on uncompressed seekable file
- biommap is tested only on Linux
- biommap is still in developpement many thing can change or be break

## Usage

### From source

In your `Cargo.toml` put
```toml
biommap = { git = "https://github.com/natir/biommap/", branch = "main" }
```

### Features

- fasta: Include fasta parser
- fastq: Include fastq parser
- vcf: Include vcf parser
- shared_state: Add shared\_state parser and `rayon` crates
- macro: WIP

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.74.
