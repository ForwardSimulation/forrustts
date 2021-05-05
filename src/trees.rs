/// Error type related to [``TreeSequence``]
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum TreesError {
    /// Raised by [``TreeSequence::new``].
    #[error("Tables not indexed.")]
    TablesNotIndexed,
}

pub struct TreeSequence<'a> {
    tables: &'a crate::TableCollection,
}

/// Result type for operations on trees and tree sequences.
pub type TreesResult<T> = Result<T, TreesError>;

impl<'a> TreeSequence<'a> {
    pub fn new(tables: &'a crate::TableCollection) -> TreesResult<Self> {
        Ok(Self { tables })
    }

    pub fn tables(&self) -> &'a crate::TableCollection {
        self.tables
    }
}

#[cfg(test)]
mod test_trees {

    #[test]
    fn test_treeseq_creation_and_table_access() {
        let mut tables = crate::TableCollection::new(100).unwrap();
        tables.add_edge(0, 1, 0, 1).unwrap();

        let ts = super::TreeSequence::new(&tables).unwrap();

        let tref = ts.tables();
        assert_eq!(tref.edges().len(), 1);
    }
}
