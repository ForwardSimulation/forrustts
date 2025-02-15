//! Core types for [`forrustts`](https://docs.rs/forrustts)

#![warn(missing_docs)]
#![warn(rustdoc::broken_intra_doc_links)]

use thiserror::Error;

mod position;
pub mod prelude;
#[cfg(feature = "rand")]
mod rand_position;
mod time;

pub use position::Position;
pub use time::Time;

/// Error type
#[derive(Error, Debug)]
pub enum Error {
    /// Invalid [`Position`]
    ///
    /// # Example
    ///
    /// ```
    /// if let Err(e) = forrustts_core::Position::try_from(-1) {
    ///     assert!(matches!(e, forrustts_core::Error::PositionError(-1)));
    /// }
    /// # else {
    /// #   panic!("should be an Err...");
    /// # }
    /// ```
    #[error("{0:?}")]
    PositionError(i64),
}
