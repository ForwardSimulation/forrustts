#![warn(missing_docs)]

//! Rust library for forward time population
//! genetic simulation with tree sequence recording.
//!
//! # Overview
//!
//! This package is a port of many ideas from
//! [fwdpp](https://github.com/molpopgen/fwdpp) from C++ to rust.
//!
//! Currently, this should be viewed as **experimental**.
//!
//! # Where to find examples
//!
//! In the `examples/` directory of the project repository.

pub use forrustts_core::newtypes::*;
pub use forrustts_simplification::*;
pub use forrustts_tables::*;
pub use forrustts_trees::*;

pub mod tskit_export;

/// Get the forrustts version number.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
