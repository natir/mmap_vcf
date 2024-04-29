//! Fastq parsing

/* std use */

/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::block;
use crate::error;

#[cfg(feature = "derive")]
#[biommap_derive::file2block(name = File2Block, block_type = memmap2::Mmap)]
fn fastq(block: &[u8]) -> error::Result<u64> {
    let mut end = block.len();

    for _ in 0..5 {
        end = block[..end]
            .rfind_byte(b'\n')
            .ok_or(error::Error::NoNewLineInBlock)?;

        if end + 1 < block.len() && block[end + 1] == b'@' {
            let prev = block[..end]
                .rfind_byte(b'\n')
                .ok_or(error::Error::NoNewLineInBlock)?;
            if block[prev + 1] == b'+' {
                let prevprev = block[..prev]
                    .rfind_byte(b'\n')
                    .ok_or(error::Error::NoNewLineInBlock)?;
                if block[prevprev + 1] == b'+' {
                    return Ok((end + 1) as u64);
                } else {
                    let prevprevprev = block[..prevprev]
                        .rfind_byte(b'\n')
                        .ok_or(error::Error::NoNewLineInBlock)?;
                    if block[prevprevprev + 1] == b'@' {
                        return Ok((prevprevprev + 1) as u64);
                    } else {
                        return Err(error::Error::NotAFastqFile);
                    }
                }
            } else {
                return Ok((end + 1) as u64);
            }
        }
    }

    Err(error::Error::NotAFastqFile)
}

/// Struct that store a fastq record
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Record<'a> {
    comment: &'a [u8],
    sequence: &'a [u8],
    plus: &'a [u8],
    quality: &'a [u8],
}

impl<'a> Record<'a> {
    /// Fastq comment, without `>`
    pub fn comment(&self) -> &'a [u8] {
        &self.comment
    }

    /// Fastq sequence
    pub fn sequence(&self) -> &'a [u8] {
        &self.sequence
    }

    /// Fastq plus line, without `+`
    pub fn plus(&self) -> &'a [u8] {
        &self.plus
    }

    /// Fastq quality
    pub fn quality(&self) -> &'a [u8] {
        &self.quality
    }
}

#[cfg(feature = "derive")]
#[biommap_derive::block2record(name = Block2Record, generic_type = DATA)]
pub fn fastq(&mut self) -> error::Result<Option<Record<'_>>> {
    if self.offset == self.block.len() {
        Ok(None)
    } else {
        let comment = &self.block.data()[self.get_line()?];
        self.offset += comment.len() as u64 + 1;

        let sequence = &self.block.data()[self.get_line()?];
        self.offset += sequence.len() as u64 + 1;

        let plus = &self.block.data()[self.get_line()?];
        self.offset += plus.len() as u64 + 1;

        let quality = &self.block.data()[self.get_line()?];
        self.offset += quality.len() as u64 + 1;

        Ok(Some(Record {
            comment,
            sequence,
            plus,
            quality,
        }))
    }
}

#[cfg(test)]
mod tests {
    /* std use */
    use std::io::Write as _;

    /* crates use */
    use biotest::Format as _;

    /* project use */

    /* local use */
    use super::*;

    #[cfg(feature = "derive")]
    #[test]
    fn default() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().sequence_len(150).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut blocks = File2Block::new(temp.path())?;

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 7866);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 7866);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 7866);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 7866);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 2736);

        let option = blocks.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn blocksize() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().sequence_len(150).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut blocks = File2Block::with_blocksize(8192 * 2, temp.path())?;

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 16074);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 16074);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 2052);

        let option = blocks.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn offset() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().sequence_len(150).build().unwrap();

        generator.create(temp.path(), &mut rng, 10).unwrap();

        let mut blocks = File2Block::with_offset(1000, temp.path())?;

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 2420);

        let option = blocks.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn blocksize_offset() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().sequence_len(150).build().unwrap();

        generator.create(temp.path(), &mut rng, 100).unwrap();

        let mut blocks = File2Block::with_blocksize_offset(4096, 24014, temp.path())?;

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 4030);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 3762);

        let option = blocks.next_block()?;
        assert!(option.is_some());
        let block = option.unwrap();
        assert_eq!(block.len(), 2394);

        let option = blocks.next_block()?;
        assert!(option.is_none());

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn records() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut rng = biotest::rand();
        let generator = biotest::Fastq::builder().sequence_len(50).build().unwrap();

        generator.create(temp.path(), &mut rng, 10).unwrap();

        let mut comments = Vec::new();
        let mut seqs = Vec::new();
        let mut pluss = Vec::new();
        let mut quals = Vec::new();

        let mut blocks = File2Block::new(temp.path())?;

        while let Ok(Some(block)) = blocks.next_block() {
            let mut reader = Block2Record::new(block);

            while let Ok(Some(record)) = reader.next_record() {
                comments.push(String::from_utf8(record.comment().to_vec()).unwrap());
                seqs.push(String::from_utf8(record.sequence().to_vec()).unwrap());
                pluss.push(String::from_utf8(record.plus().to_vec()).unwrap());
                quals.push(String::from_utf8(record.quality().to_vec()).unwrap());
            }
        }

        assert_eq!(
            comments,
            vec![
                "@GSWNPZYBHL atbbutlfemxuzgaghmwn".to_string(),
                "@LHABRFOIPY kuougqcoherrilpylxga".to_string(),
                "@PHYALDFVGY kzhdadxmhhvhzlvxulto".to_string(),
                "@UWWIFBAZNV hmqfbhzkizjnwzyvreoi".to_string(),
                "@ETFPQPBCGR uiwejacxekdfeoxhehyv".to_string(),
                "@OXRTYXGRDF nagvxtzlfvryzktzbvur".to_string(),
                "@WWMYBRMIVZ utdjngxdvdwshbpjuvvr".to_string(),
                "@JIGLPJJGED cphroskdsgmlwqixchga".to_string(),
                "@DEPZJRDOQP mdadqdzktocmqeqfubqc".to_string(),
                "@DKYULMZDLL qgqrqlfpzyvfawnwhibf".to_string(),
            ]
        );
        assert_eq!(
            seqs,
            vec![
                "gccAcggtAatGcTtgtaCgcAGgAtaTcgAAtTaTaGaTggttGCtCat".to_string(),
                "ActTgCGAcaAgaAAtaTCCcAgagggaCcttCcGcTTGcgAACcTtCtt".to_string(),
                "GCatAGGACAAaacTaTTagagGtatAGCcTatTtaaaaCGgcttGGTtg".to_string(),
                "CcCGtCtATgTTgTATcaTTCGaCCttcAaGCGCAatgaTGAtaatcaCt".to_string(),
                "caacCATAtGtGgtAtacAagtTGgAtgcGtTCtctTgctTtcGggATtc".to_string(),
                "AttaTtGgCACgCcgcCgATtcGCaTatTGGGCTacgtgACCGttTCAtt".to_string(),
                "cTGAAccCaaAcacagCATCTaTCgGcgcaGCaCaTATTacCGaTtgttC".to_string(),
                "CatTagaCtTTtccCCcAgagtctAGCCtCTgATTtTGCcGcGgCgTcGc".to_string(),
                "gGtAGgATcAaGactACcGCAaatGtTgCtaccTaCGCgGGaacgaCcAt".to_string(),
                "TCatCGgcAgaaagtTACacgcggaCcTAccaaGgGCAAGtGcccCgaGc".to_string(),
            ]
        );
        assert_eq!(
            pluss,
            vec![
                "+Wognt".to_string(),
                "+gGXi[".to_string(),
                "+T_IlS".to_string(),
                "+eGSa[".to_string(),
                "+YcBuG".to_string(),
                "+P_hBj".to_string(),
                "+jXNUP".to_string(),
                "+R`mig".to_string(),
                "+SQrSl".to_string(),
                "+NS^fl".to_string(),
            ]
        );
        assert_eq!(
            quals,
            vec![
                "/&411+!)AF,E;7.8.3GF2%\"%:4%#<399BE%$8900(08#,.;&2*".to_string(),
                ">B4-?CG@!GI!(?\"1'%))7&08<27F?3AA$E(/@A#FBBF<')G+%2".to_string(),
                "?\"))&?,%6;GG@D-.4G5+%$8H>&%@72B;''I:>+.*\"G.$8E&./=".to_string(),
                "!20/(B?(5\"/FC.D7*-I,%#G%;4\"7->$D@6=,C2\"')7G0FA-+I'".to_string(),
                "9?A)>5,>E6G9<E+4@)D+\")9@>;B3#.FG;:.#21\"I378/0G=377".to_string(),
                "17?)A$AA#C\"DI%H3F/9;'H1G93#'+4<AF/$*B,52-@!C=\"G<63".to_string(),
                "8,%$'\"6.>:#-'G1BH11)&97>HI0D#(96H.0D0G1'(FDG::-;4\"".to_string(),
                "7C\"&C?\"=A46?D9I1=<1345'/-;%F6363*70@$1;:80>497,@1>".to_string(),
                "$.G-\"4?(-H&%,7<>;EE4F%00#5H2B45;2*/=I:;=36H5$%C(08".to_string(),
                "4C+*#).<?3?1?3.??.=E:&*--;G2.?AIC('D$81H@5A;=$*D<3".to_string(),
            ]
        );

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn quality_is_shit() -> error::Result<()> {
        let data = b"@1\nAA\n+1\n!!\n@2\nTT\n+2\n!!";
        assert_eq!(File2Block::correct_block_size(data)?, 12);

        let data = b"@1\nAA\n+1\n!!\n@2\nTT\n+2\n+!\n@3";
        assert_eq!(File2Block::correct_block_size(data)?, 24);

        let data = b"@1\nAA\n+1\n!!\n@2\nTT\n+2\n@!";
        assert_eq!(File2Block::correct_block_size(data)?, 12);

        Ok(())
    }

    #[cfg(feature = "derive")]
    #[test]
    fn not_a_fastq() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;

        let mut data = b"@0
TTAGATTATAGTACGG
ATTATAT
+1
AGTTATCGTGTACCTC
+1
+CW?:KL~15\\E|MN
GTCCCTCAATCCG
+2
"
        .to_vec();

        {
            let mut file = std::io::BufWriter::new(std::fs::File::create(&temp.path())?);
            file.write_all(&data)?;
        }

        let mut blocks = File2Block::with_blocksize(82, temp.path())?;

        println!("{:?} {}", String::from_utf8(data.clone()), data.len());
        let result = blocks.next_block();
        println!("{:?}", result);
        assert!(result.is_err());

        data.extend(
            b"+FAILLED FILE
+3
+TTGGGCATGAGGTTCA
@3ueauie
+~vGLKg+n!*iJ\\K
@iuiea
",
        );

        {
            let mut file = std::io::BufWriter::new(std::fs::File::create(&temp.path())?);
            file.write_all(&data)?;
        }

        let mut blocks = File2Block::with_blocksize(82, temp.path())?;

        assert!(blocks.next_block().is_err());

        let mut blocks = File2Block::with_blocksize(82, temp.path())?;
        assert!(blocks.next().is_some());
        assert!(blocks.next().unwrap().is_err());

        Ok(())
    }
}
