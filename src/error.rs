use crate::nested_forward_list::NestedForwardListError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ForrusttsError {
    #[error("{value:?}")]
    SimplificationError { value: String },
    #[error("{value:?}")]
    ListError {
        #[from]
        value: NestedForwardListError,
    },
    #[error("{value:?}")]
    TablesError {
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
            Err(NestedForwardListError::InvalidKey)
        }
    }

    fn return_simplification_error() -> Result<(), ForrusttsError> {
        match return_nested_forward_list_error(false) {
            Ok(_) => Ok(()),
            Err(e) => match e {
                NestedForwardListError::InvalidKey => Err(ForrusttsError::ListError { value: e }),
                NestedForwardListError::NullTail => Err(ForrusttsError::ListError { value: e }),
                NestedForwardListError::KeyOutOfRange => {
                    Err(ForrusttsError::ListError { value: e })
                }
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
                    NestedForwardListError::InvalidKey => {
                        assert_eq!(value.to_string(), "Invalid key")
                    }
                    NestedForwardListError::NullTail => panic!(),
                    NestedForwardListError::KeyOutOfRange => panic!(),
                },
                #[allow(unreachable_patterns)]
                _ => panic!(),
            },
        };
    }
}
