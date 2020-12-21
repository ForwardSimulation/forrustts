mod edge_buffer;
pub mod nested_forward_list;
mod segment;
mod simplification_buffers;
mod simplification_flags;
mod simplification_logic;
mod simplification_output;
mod simplify_from_edge_buffer;
mod simplify_tables;
mod tables;
mod tsdef;

pub use edge_buffer::EdgeBuffer;
pub use segment::Segment;
pub use simplification_buffers::SimplificationBuffers;
pub use simplification_flags::SimplificationFlags;
pub use simplification_output::SimplificationOutput;
pub use simplify_from_edge_buffer::simplify_from_edge_buffer;
pub use simplify_tables::{simplify_tables, simplify_tables_with_buffers};
pub use tables::*;
pub use tsdef::*;

pub mod tskit;
pub mod wright_fisher;

/// Get the tskit_rust version number.
pub fn version() -> &'static str {
    return env!("CARGO_PKG_VERSION");
}

// These are testing modules
mod test_simplify_tables;
