mod macros;

pub mod newtypes;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("range error: {}", *.0)]
    ConversionError(String),
}
