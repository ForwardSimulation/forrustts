///! Error handling
use crate::nested_forward_list::NestedForwardListError;
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

    /// A redirection of a [``crate::nested_forward_list::NestedForwardListError``].
    #[error("{value:?}")]
    ListError {
        /// The redirected error
        #[from]
        value: NestedForwardListError,
    },
    /// A redirection of a [``crate::TablesError``]
    #[error("{value:?}")]
    TablesError {
        /// The redirected error
        #[from]
        value: crate::TablesError,
    },
}

#[cfg(test)]
mod test {

    use super::*;

    fn return_nested_forward_list_error(f: bool) -> Result<(), NestedForwardListError> {
        if f {
            Ok(())
        } else {
            Err(NestedForwardListError::InvalidIndex)
        }
    }

    fn return_simplification_error() -> Result<(), ForrusttsError> {
        match return_nested_forward_list_error(false) {
            Ok(_) => Ok(()),
            Err(e) => match e {
                NestedForwardListError::InvalidIndex => Err(ForrusttsError::ListError { value: e }),
                NestedForwardListError::NullTail => Err(ForrusttsError::ListError { value: e }),
                #[allow(unreachable_patterns)]
                _ => panic!(),
            },
        }
    }

    #[test]
    fn test_nested_forward_list_error_propagation() {
        match return_simplification_error() {
            Ok(_) => panic!(),
            Err(e) => match e {
                ForrusttsError::ListError { value } => match value {
                    NestedForwardListError::InvalidIndex => {
                        assert_eq!(value.to_string(), "Invalid index")
                    }
                    NestedForwardListError::NullTail => panic!(),
                },
                #[allow(unreachable_patterns)]
                _ => panic!(),
            },
        };
    }
}
