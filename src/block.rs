//! Define a Block of file

/* std use */

/* crate use */

/* project use */

/* mod declaration */
#[cfg(derive)]
pub mod derive;

/// Manage a block of record in DATA indexable struct
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Block<DATA> {
    data: DATA,
    length: u64,
}

impl<DATA> Block<DATA> {
    /// Create a new Block
    pub fn new(data: DATA, length: u64) -> Self {
        Self { data, length }
    }

    /// Get length of block
    pub fn len(&self) -> u64 {
        self.length
    }

    /// Return true if the block is empty
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }
}

impl<DATA> Block<DATA>
where
    DATA: core::convert::AsRef<[u8]>,
{
    /// Acces to data owned by block
    pub fn data(&self) -> &[u8] {
        &self.data.as_ref()[..self.length as usize]
    }

    /// Return remain part of block
    pub fn remain(&self) -> &[u8] {
        &self.data.as_ref()[self.length as usize..]
    }
}

#[cfg(test)]
mod tests {
    /* std use */
    use core::ops::Deref;

    /* crate use */
    use biotest::Format as _;

    /* project use */
    use crate::error;

    /* local use */
    use super::*;

    #[test]
    fn ram() -> error::Result<()> {
        let mut rng = biotest::rand();
        let mut data = Vec::new();
        let generator = biotest::Fasta::builder().build().unwrap();

        generator.records(&mut data, &mut rng, 5).unwrap();

        let block = Block::new(data.clone(), 300);

        assert_eq!(block.data(), &data[..300]);
        assert_eq!(block.remain(), &data[300..]);
        assert_eq!(block.remain().len(), 620);
        assert!(!block.is_empty());

        Ok(())
    }

    #[test]
    fn mmap() -> error::Result<()> {
        let temp = tempfile::NamedTempFile::new()?;
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().build().unwrap();

        generator.create(temp.path(), &mut rng, 5).unwrap();

        let data = unsafe {
            memmap2::MmapOptions::new()
                .offset(0)
                .len(500)
                .map(temp.as_file())?
        };

        let block = Block::new(data.deref(), 300);

        assert_eq!(block.data(), &data[..300]);
        assert_eq!(block.remain(), &data[300..]);
        assert_eq!(block.remain().len(), 200);
        assert!(!block.is_empty());

        Ok(())
    }

    #[test]
    fn empty() -> error::Result<()> {
        let data = vec![0; 0];

        let block = Block::new(data.clone(), 0);

        assert_eq!(block.data(), &data[..0]);
        assert_eq!(block.remain(), &data[0..]);
        assert_eq!(block.remain().len(), 0);
        assert!(block.is_empty());

        Ok(())
    }
}
