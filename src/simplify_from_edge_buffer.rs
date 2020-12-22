use crate::simplification_logic;
use crate::tables::*;
use crate::EdgeBuffer;
use crate::ForrusttsError;
use crate::Segment;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;
use crate::{IdType, Position, Time};

struct ParentLocation {
    parent: IdType,
    start: usize,
    stop: usize,
}

// TODO: validate input and return errors.
impl ParentLocation {
    fn new(parent: IdType, start: usize, stop: usize) -> Self {
        ParentLocation {
            parent,
            start,
            stop,
        }
    }
}

fn find_pre_existing_edges(
    tables: &TableCollection,
    alive_at_last_simplification: &[IdType],
    edge_buffer: &EdgeBuffer,
) -> Result<Vec<ParentLocation>, ForrusttsError> {
    let mut alive_with_new_edges: Vec<i32> = vec![];

    for a in alive_at_last_simplification {
        if edge_buffer.head(*a)? != EdgeBuffer::null() {
            alive_with_new_edges.push(*a);
        }
    }
    if alive_with_new_edges.is_empty() {
        return Ok(vec![]);
    }

    let mut starts = vec![usize::MAX; tables.num_nodes()];
    let mut stops = vec![usize::MAX; tables.num_nodes()];

    for (i, e) in tables.enumerate_edges() {
        if starts[e.parent as usize] == usize::MAX {
            starts[e.parent as usize] = i;
            stops[e.parent as usize] = i + 1;
        } else {
            stops[e.parent as usize] = i + 1;
        }
    }

    let mut rv = vec![];
    for a in alive_with_new_edges {
        rv.push(ParentLocation::new(
            a,
            starts[a as usize],
            stops[a as usize],
        ));
    }

    rv.sort_by(|a, b| {
        let ta = tables.nodes_[a.parent as usize].time;
        let tb = tables.nodes_[b.parent as usize].time;
        if ta == tb {
            if a.start == b.start {
                return a.parent.cmp(&b.parent);
            }
            return a.start.cmp(&b.start);
        }
        ta.cmp(&tb).reverse()
    });

    // TODO: this could eventually be called in a debug_assert
    if !rv.is_empty() {
        for i in 1..rv.len() {
            let t0 = tables.nodes_[rv[i - 1].parent as usize].time;
            let t1 = tables.nodes_[rv[i].parent as usize].time;
            if t0 < t1 {
                return Err(ForrusttsError::SimplificationError {
                    value: "existing edges not properly sorted by time".to_string(),
                });
            }
        }
    }
    Ok(rv)
}

fn queue_children(
    child: IdType,
    left: Position,
    right: Position,
    ancestry: &mut simplification_logic::AncestryList,
    overlapper: &mut simplification_logic::SegmentOverlapper,
) -> Result<(), ForrusttsError> {
    Ok(ancestry.for_each(child, |seg: &Segment| {
        if seg.right > left && right > seg.left {
            overlapper.enqueue(
                std::cmp::max(seg.left, left),
                std::cmp::min(seg.right, right),
                seg.node,
            );
        }
        true
    })?)
}

fn process_births_from_buffer(
    head: IdType,
    edge_buffer: &EdgeBuffer,
    state: &mut SimplificationBuffers,
) -> Result<(), ForrusttsError> {
    // Have to take references here to
    // make the borrow checker happy.
    let a = &mut state.ancestry;
    let o = &mut state.overlapper;
    Ok(edge_buffer.for_each(head, |seg: &Segment| {
        queue_children(seg.node, seg.left, seg.right, a, o).unwrap();
        true
    })?)
}

pub fn simplify_from_edge_buffer(
    samples: &[IdType],
    alive_at_last_simplification: &[IdType],
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    edge_buffer: &mut EdgeBuffer,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), ForrusttsError> {
    // FIXME: validate alive_at_last_simplification

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

    simplification_logic::setup_idmap(&tables.nodes_, &mut output.idmap);

    state.clear();
    state.ancestry.reset(tables.num_nodes());

    simplification_logic::record_sample_nodes(
        &samples,
        &tables,
        &mut state.new_nodes,
        &mut state.ancestry,
        &mut output.idmap,
    )?;

    // Process all edges since the last simplification.
    let mut max_time = Time::MIN;
    for n in alive_at_last_simplification {
        max_time = std::cmp::max(max_time, tables.node(*n).time);
    }
    for (i, _) in edge_buffer.head_itr().rev().enumerate() {
        let head = (edge_buffer.len() - i - 1) as i32;
        let ptime = tables.node(head).time;
        if ptime > max_time
        // Then this is a parent who is:
        // 1. Born since the last simplification.
        // 2. Left offspring
        {
            state.overlapper.clear_queue();
            process_births_from_buffer(head, edge_buffer, state)?;
            state.overlapper.finalize_queue(tables.get_length());
            simplification_logic::merge_ancestors(
                &tables.nodes(),
                tables.get_length(),
                head,
                state,
                &mut output.idmap,
            )?;
        } else if ptime <= max_time {
            break;
        }
    }

    let existing_edges =
        find_pre_existing_edges(&tables, &alive_at_last_simplification, &edge_buffer)?;

    let mut edge_i = 0;
    let num_edges = tables.num_edges();

    for ex in existing_edges {
        while edge_i < num_edges
            && tables.nodes_[tables.edges_[edge_i].parent as usize].time
                > tables.nodes_[ex.parent as usize].time
        {
            let u = tables.edges_[edge_i].parent;
            edge_i = simplification_logic::find_parent_child_segment_overlap(
                &tables.edges_,
                edge_i,
                num_edges,
                tables.get_length(),
                u,
                &mut state.ancestry,
                &mut state.overlapper,
            )?;
            simplification_logic::merge_ancestors(
                &tables.nodes_,
                tables.get_length(),
                u,
                state,
                &mut output.idmap,
            )?;
        }
        if ex.start != usize::MAX {
            while (edge_i as usize) < ex.start
                && tables.nodes_[tables.edges_[edge_i].parent as usize].time
                    >= tables.nodes_[ex.parent as usize].time
            {
                let u = tables.edges_[edge_i].parent;
                edge_i = simplification_logic::find_parent_child_segment_overlap(
                    &tables.edges_,
                    edge_i,
                    num_edges,
                    tables.get_length(),
                    u,
                    &mut state.ancestry,
                    &mut state.overlapper,
                )?;
                simplification_logic::merge_ancestors(
                    &tables.nodes_,
                    tables.get_length(),
                    u,
                    state,
                    &mut output.idmap,
                )?;
            }
        }
        // now, handle ex.parent
        state.overlapper.clear_queue();
        if ex.start != usize::MAX {
            while edge_i < ex.stop {
                // TODO: a debug assert or regular assert?
                if tables.edges_[edge_i].parent != ex.parent {
                    return Err(ForrusttsError::SimplificationError {
                        value: "Unexpected parent node".to_string(),
                    });
                }
                let a = &mut state.ancestry;
                let o = &mut state.overlapper;
                queue_children(
                    tables.edges_[edge_i].child,
                    tables.edges_[edge_i].left,
                    tables.edges_[edge_i].right,
                    a,
                    o,
                )?;
                edge_i += 1;
            }
            if edge_i < num_edges && tables.edges_[edge_i].parent == ex.parent {
                return Err(ForrusttsError::SimplificationError {
                    value: "error traversing pre-existing edges for parent".to_string(),
                });
            }
        }
        process_births_from_buffer(ex.parent, edge_buffer, state)?;
        state.overlapper.finalize_queue(tables.get_length());
        simplification_logic::merge_ancestors(
            &tables.nodes_,
            tables.get_length(),
            ex.parent,
            state,
            &mut output.idmap,
        )?;
    }

    // Handle remaining edges.
    while edge_i < num_edges {
        let u = tables.edges_[edge_i].parent;
        edge_i = simplification_logic::find_parent_child_segment_overlap(
            &tables.edges_,
            edge_i,
            num_edges,
            tables.get_length(),
            u,
            &mut state.ancestry,
            &mut state.overlapper,
        )?;

        simplification_logic::merge_ancestors(
            &tables.nodes_,
            tables.get_length(),
            u,
            state,
            &mut output.idmap,
        )?;
    }

    std::mem::swap(&mut tables.edges_, &mut state.new_edges);
    std::mem::swap(&mut tables.nodes_, &mut state.new_nodes);
    edge_buffer.reset(tables.num_nodes());

    Ok(())
}

#[cfg(test)]
mod test {

    use super::*;

    // This shows that the closure error gets propagated
    // as the result type.
    #[test]
    fn test_process_births_from_buffer_closure_error() {
        let b = EdgeBuffer::new();
        let mut s = SimplificationBuffers::new();
        assert!(process_births_from_buffer(-1, &b, &mut s)
            .map_or_else(|_: ForrusttsError| true, |_| false));
    }
}
