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

// FIXME: we are confusing necessary exports with
// what should be in a prelude?
pub use forrustts_core::newtypes::*;
pub use forrustts_tables_trees::prelude::*;
pub use forrustts_tskit::export_tables;
pub use forrustts_tskit::TableCollectionExportFlags;
