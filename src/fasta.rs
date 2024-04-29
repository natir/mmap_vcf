//! Fasta parsing

/* std use */

/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::block;
use crate::error;

#[cfg(feature = "derive")]
#[biommap_derive::file2block(name = File2Block, block_type = memmap2::Mmap)]
fn fasta(block: &[u8]) -> error::Result<u64> {
    let mut end = block.len();

    for _ in 0..2 {
        end = block[..end]
            .rfind_byte(b'\n')
            .ok_or(error::Error::NoNewLineInBlock)?;

        if end + 1 < block.len() && block[end + 1] == b'>' {
            return Ok((end + 1) as u64);
        }
    }

    Err(error::Error::NotAFastaFile)
}

/// Struct that store a fasta record
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Record<'a> {
    comment: &'a [u8],
    sequence: &'a [u8],
}

impl<'a> Record<'a> {
    /// Fasta comment without `>`
    pub fn comment(&self) -> &'a [u8] {
        &self.comment
    }

    /// Fasta sequence
    pub fn sequence(&self) -> &'a [u8] {
        &self.sequence
    }
}

#[cfg(feature = "derive")]
#[biommap_derive::block2record(name = Block2Record, generic_type = DATA)]
pub fn fasta(&mut self) -> error::Result<Option<Record<'_>>> {
    if self.offset == self.block.len() {
        Ok(None)
    } else {
        let comment = &self.block.data()[self.get_line()?];
        self.offset += comment.len() as u64 + 1;

        let sequence = &self.block.data()[self.get_line()?];
        self.offset += sequence.len() as u64 + 1;

        Ok(Some(Record { comment, sequence }))
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* crates use */
    use biotest::Format as _;

    /* project use */
    use crate::error;

    /* local use */
    use super::*;

    #[cfg(feature = "derive")]
    #[test]
    fn default() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(50).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut producer = File2Block::new(temp.path())?;

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 8148);

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 252);

        let option = producer.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn blocksize() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(50).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut producer = File2Block::with_blocksize(8192 * 2, temp.path())?;

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 8400);

        let option = producer.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn offset() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(50).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut producer = File2Block::with_offset(8050, temp.path())?;

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 350);

        let option = producer.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn blocksize_offset() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(100).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut producer = File2Block::with_blocksize_offset(4096, 8050, temp.path())?;

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 4010);

        let option = producer.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 1340);

        let option = producer.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn records() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(50).build().unwrap();

        generator.create(temp.path(), &mut rng, 10).unwrap();

        let mut comments = Vec::new();
        let mut seqs = Vec::new();

        let mut producer = File2Block::new(temp.path())?;

        while let Ok(Some(block)) = producer.next_block() {
            let mut reader = Block2Record::new(block);

            while let Ok(Some(record)) = reader.next_record() {
                comments.push(String::from_utf8(record.comment().to_vec()).unwrap());
                seqs.push(String::from_utf8(record.sequence().to_vec()).unwrap());
            }
        }

        assert_eq!(
            comments,
            vec![
                ">GSWNPZYBHL atbbutlfemxuzgaghmwn".to_string(),
                ">RCUDMHKKGS ajefuxhqoiwnooilwywl".to_string(),
                ">RVOLOAYNLY acavbgerslbixoxxodry".to_string(),
                ">CLCUXEJUUO bxvvxhrfygckrphyaldf".to_string(),
                ">UNLMOTCONV zzkwrudmkpkjusxndtiw".to_string(),
                ">IZJNWZYVRE oizferdlsseuahsvxvjh".to_string(),
                ">MOHSCWTGTN hvlfyxahfdjoyxuahmga".to_string(),
                ">MELUFGTSRI ugyirugryxamshjpzprp".to_string(),
                ">AGVXTZLFVR yzktzbvurjcfibwtjutf".to_string(),
                ">ASKDOTFRUC uubdjpvcftawzzlxspaf".to_string()
            ]
        );
        assert_eq!(
            seqs,
            vec![
                "gccAcggtAatGcTtgtaCgcAGgAtaTcgAAtTaTaGaTggttGCtCat".to_string(),
                "AGacAtgCtGCAAtTacCGtTAAcaGGtatTCaTCctcTGgAActTgCGA".to_string(),
                "ttCcGcTTGcgAACcTtCttAacGtTtAtGTgACAGCCaCGctGagattT".to_string(),
                "TGTCCACgTTTGagtGaGCatAGGACAAaacTaTTagagGtatAGCcTat".to_string(),
                "ACTacgtCTaTgTCAGgCtaGTtcCCTcgcTgAgGgAtCAAatTCTATTG".to_string(),
                "ATcaTTCGaCCttcAaGCGCAatgaTGAtaatcaCtGcTAGCCAgaTTgc".to_string(),
                "CCtcTctCAtgCGCagTCTcaacCATAtGtGgtAtacAagtTGgAtgcGt".to_string(),
                "AgtaTgacgtCCTAtActaGAggcAAGGACGaATctgCaaatgctgTcCa".to_string(),
                "aTtGgCACgCcgcCgATtcGCaTatTGGGCTacgtgACCGttTCAttTac".to_string(),
                "GgACTctgTGTtaAGCAgcagAcGttCagTgCTAtccTGAAccCaaAcac".to_string(),
            ]
        );

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn not_fasta() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut producer = File2Block::with_blocksize(200, temp.path())?;

        let result = producer.next_block();
        assert_matches::assert_matches!(result, Err(error::Error::NotAFastaFile));

        Ok(())
    }
}
