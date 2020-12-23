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
//! # Differences from [tskit](https://tskit.readthedocs.io)
//!
//! The tables model differs from `tskit` in some important ways:
//!
//! 1. Time moves from the past to the present.
//!    Thus, child nodes have time values *greater than*
//!    those of their parents.
//! 2. The data layout is "array of structures" while
//!    `tskit` is a "structure of arrays".
//! 3. Metadata is not part of the tables.
//!    Often, what one thinks of as metadata is data
//!    used during the simulation.  Thus, it is not
//!    part of a [``TableCollection``] and is something
//!    that would only be useful to write when transfering
//!    final results ot a [``tskit_rust::TableCollection``].
//! 4. Mutation table data are different. See [``MutationRecord``].
//!
//! # Where to find examples
//!
//! The [repository](https://github.com/molpopgen/forrustts)
//! contains a "sub-crate" of examples, called `forrustts_examples`.
//!
//! The source code of the [``wright_fisher``] module contains a full
//! example of Wright-Fisher simulation with recombination and overlapping.
//! generations.

// NOTE: uncomment the next line in order to find
// stuff that needs documenting:
// #![warn(missing_docs)]

mod edge_buffer;
mod error;
pub mod nested_forward_list;
mod segment;
mod simplification_buffers;
mod simplification_common;
mod simplification_flags;
mod simplification_logic;
mod simplification_output;
mod simplify_from_edge_buffer;
mod simplify_tables;
mod tables;
mod tsdef;

pub use edge_buffer::EdgeBuffer;
pub use error::ForrusttsError;
pub use segment::Segment;
pub use simplification_buffers::SimplificationBuffers;
pub use simplification_flags::SimplificationFlags;
pub use simplification_output::SimplificationOutput;
pub use simplify_from_edge_buffer::simplify_from_edge_buffer;
pub use simplify_tables::{simplify_tables, simplify_tables_without_state};
pub use tables::*;
pub use tsdef::*;

pub mod tskit;
pub mod wright_fisher;

/// Get the forrustts version number.
pub fn version() -> &'static str {
    return env!("CARGO_PKG_VERSION");
}

// These are testing modules
mod test_simplify_tables;
