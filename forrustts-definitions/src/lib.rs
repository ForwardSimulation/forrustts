//! A light crate defining types and constants.

/// Integer type used to refer to a genomic coordinate/position
pub type Position = i64;

/// Time is recorded as a discrete, signed integer.
pub type Time = i64;

/// Integer type used to refer to objects by their "id".
pub type IdType = i32;

/// Equals -1 (minus one).
/// Primary use is to indicate a null "id" value.
pub static NULL_ID: IdType = -1;
