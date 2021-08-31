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
//!    final results to a `tskit::TableCollection`
//!    (see [`tskit`](https://crates.io/crates/tskit)).
//! 4. Mutation table data are different. See [``MutationRecord``].
//! 5. Genomic locations are integers (see [``Position``]).
//!    In `tskit`, both are C `double`, the equivalent of [``f64``].
//!
//! # Where to find examples
//!
//! In the `examples/` directory of the project repository.

// NOTE: uncomment the next line in order to find
// stuff that needs documenting:
// #![warn(missing_docs)]

mod macros;

pub mod nested_forward_list;
mod newtypes;
mod segment;
mod simplification;
mod tables;
mod traits;
mod trees;

pub use newtypes::*;
use segment::Segment;
pub use simplification::simplify_from_edge_buffer;
pub use simplification::EdgeBuffer;
pub use simplification::SamplesInfo;
pub use simplification::SimplificationBuffers;
pub use simplification::SimplificationError;
pub use simplification::SimplificationFlags;
pub use simplification::SimplificationOutput;
pub use simplification::{simplify_tables, simplify_tables_without_state};
pub use tables::*;
pub use traits::*;
pub use trees::*;
pub mod prelude;

/// Get the forrustts version number.
pub fn version() -> &'static str {
    return env!("CARGO_PKG_VERSION");
}
