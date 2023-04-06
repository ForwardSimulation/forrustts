use thiserror::Error;

mod position;
pub mod prelude;
#[cfg(feature = "rand")]
mod rand_position;
mod time;

pub use position::Position;
pub use time::Time;

#[derive(Error, Debug)]
pub enum Error {
    #[error("{0:?}")]
    PositionError(i64),
}
