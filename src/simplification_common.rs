/// Common functions to reuse in various "simplify tables"
/// functions
use crate::simplification_logic;
use crate::ForrusttsError;
use crate::SamplesInfo;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;
use crate::{IdType, NULL_ID};
use crate::{Node, TableCollection};

fn setup_idmap(nodes: &[Node], idmap: &mut Vec<IdType>) {
    idmap.resize(nodes.len(), NULL_ID);
    idmap.iter_mut().for_each(|x| *x = NULL_ID);
}

pub fn setup_simplification(
    samples: &SamplesInfo,
    tables: &TableCollection,
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    output: &mut SimplificationOutput,
) -> Result<(), ForrusttsError> {
    if !tables.sites_.is_empty() || !tables.mutations_.is_empty() {
        return Err(ForrusttsError::SimplificationError {
            value: "mutation simplification not yet implemented".to_string(),
        });
    }

    if flags.bits() != 0 {
        return Err(ForrusttsError::SimplificationError {
            value: "SimplificationFlags must be zero".to_string(),
        });
    }

    setup_idmap(&tables.nodes_, &mut output.idmap);

    state.clear();
    state.ancestry.reset(tables.num_nodes());

    simplification_logic::record_sample_nodes(
        &samples.samples,
        &tables,
        &mut state.new_nodes,
        &mut state.ancestry,
        &mut output.idmap,
    )?;

    Ok(())
}

pub fn process_parent(
    u: IdType,
    (edge_index, num_edges): (usize, usize),
    tables: &TableCollection,
    state: &mut SimplificationBuffers,
    output: &mut SimplificationOutput,
) -> Result<usize, ForrusttsError> {
    let edge_i = simplification_logic::find_parent_child_segment_overlap(
        &tables.edges_,
        edge_index,
        num_edges,
        tables.genome_length(),
        u,
        &mut state.ancestry,
        &mut state.overlapper,
    )?;

    simplification_logic::merge_ancestors(
        &tables.nodes_,
        tables.genome_length(),
        u,
        state,
        &mut output.idmap,
    )?;
    Ok(edge_i)
}
