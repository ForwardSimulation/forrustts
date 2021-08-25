use thiserror::Error;

/// Primary error type.
///
/// Some members of this enum implement ``From``
/// in order to redirect other error types.
#[derive(Error, Debug, PartialEq)]
pub enum ForrusttsError {
    /// An error that occurs during simplification.
    #[error("{value:?}")]
    SimplificationError {
        /// The error message
        value: String,
    },

    /// A redirection of a [``forrustts_tables_trees::nested_forward_list::NestedForwardListError``].
    #[error("{value:?}")]
    ListError {
        /// The redirected error
        #[from]
        value: forrustts_tables_trees::nested_forward_list::NestedForwardListError,
    },
    /// A redirection of a [``crate::TablesError``]
    #[error("{value:?}")]
    TablesError {
        /// The redirected error
        #[from]
        value: forrustts_tables_trees::TablesError,
    },
    /// A redirection of a [``crate::TreesError``]
    #[error("{value:?}")]
    TreesError {
        /// The redirected error
        #[from]
        value: forrustts_tables_trees::TreesError,
    },
}
