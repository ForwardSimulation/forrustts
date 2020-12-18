use crate::nested_forward_list::NestedForwardList;
use crate::segment::Segment;
use crate::simplification_buffers::SimplificationBuffers;
use crate::tables::*;
use crate::tsdef::{IdType, Position, NULL_ID};

pub struct SegmentOverlapper {
    segment_queue: Vec<Segment>,
    overlapping: Vec<Segment>,
    left: Position,
    right: Position,
    qbeg: usize,
    qend: usize,
    obeg: usize,
    oend: usize,
}

impl SegmentOverlapper {
    fn set_partition(&mut self) -> Position {
        let mut tright = Position::MAX;
        let mut b: usize = 0;

        for i in 0..self.oend {
            if self.overlapping[i].right > self.left {
                self.overlapping[b] = self.overlapping[i];
                tright = std::cmp::min(tright, self.overlapping[b].right);
                b += 1;
            }
        }

        self.oend = b;

        tright
    }

    fn num_overlaps(&self) -> usize {
        assert!(
            self.oend - self.obeg <= self.overlapping.len(),
            format!(
                "overlap details = {} {} {}",
                self.oend,
                self.obeg,
                self.overlapping.len()
            )
        );
        self.oend - self.obeg
    }

    // Public interface below

    pub const fn new() -> SegmentOverlapper {
        SegmentOverlapper {
            segment_queue: vec![],
            overlapping: vec![],
            left: 0,
            right: Position::MAX,
            qbeg: std::usize::MAX,
            qend: std::usize::MAX,
            obeg: std::usize::MAX,
            oend: std::usize::MAX,
        }
    }

    pub fn init(&mut self) {
        self.qbeg = 0;
        self.qend = self.segment_queue.len() - 1;
        assert!(self.qend < self.segment_queue.len());
        self.obeg = 0;
        self.oend = 0;
        self.overlapping.clear();
    }

    pub fn enqueue(&mut self, left: Position, right: Position, node: IdType) {
        self.segment_queue.push(Segment { left, right, node });
    }

    pub fn finalize_queue(&mut self, maxlen: Position) {
        self.segment_queue.sort_by(|a, b| a.left.cmp(&b.left));
        self.segment_queue.push(Segment {
            left: maxlen,
            right: maxlen + 1,
            node: NULL_ID,
        });
    }

    pub fn advance(&mut self) -> bool {
        let mut rv = false;

        if self.qbeg < self.qend {
            self.left = self.right;
            let mut tright = self.set_partition();
            if self.num_overlaps() == 0 {
                self.left = self.segment_queue[self.qbeg].left;
            }
            while self.qbeg < self.qend && self.segment_queue[self.qbeg].left == self.left {
                tright = std::cmp::min(tright, self.segment_queue[self.qbeg].right);
                // NOTE: I wonder how efficient this is vs C++?
                self.overlapping
                    .insert(self.oend, self.segment_queue[self.qbeg]);
                self.oend += 1;
                self.qbeg += 1;
            }
            self.right = std::cmp::min(self.segment_queue[self.qbeg].left, tright);
            rv = true;
        } else {
            self.left = self.right;
            self.right = Position::MAX;
            let tright = self.set_partition();
            if self.num_overlaps() > 0 {
                self.right = tright;
                rv = true
            }
        }

        rv
    }

    pub fn get_left(&self) -> Position {
        self.left
    }

    pub fn get_right(&self) -> Position {
        self.right
    }

    pub fn clear_queue(&mut self) {
        self.segment_queue.clear();
    }

    pub fn overlap(&self, i: usize) -> &Segment {
        &self.overlapping[i]
    }
}

pub type AncestryList = NestedForwardList<Segment>;

// FIXME: another panic! room
pub fn find_parent_child_segment_overlap(
    edges: &[Edge],
    edge_index: usize,
    num_edges: usize,
    maxlen: Position,
    u: IdType,
    ancestry: &mut AncestryList,
    overlapper: &mut SegmentOverlapper,
) -> usize {
    overlapper.clear_queue();

    let mut i = edge_index;

    while i < num_edges && edges[i].parent == u {
        let edge = &edges[i];
        ancestry
            .consume(edges[i].child, |seg: &Segment| {
                if seg.right > edge.left && edge.right > seg.left {
                    overlapper.enqueue(
                        std::cmp::max(seg.left, edge.left),
                        std::cmp::min(seg.right, edge.right),
                        seg.node,
                    );
                }
                true
            })
            .unwrap();
        i += 1;
    }
    overlapper.finalize_queue(maxlen);
    i
}

pub fn setup_idmap(nodes: &[Node]) -> Vec<IdType> {
    return vec![NULL_ID; nodes.len()];
}

// FIXME: this will panic! if we get error from AncestryList
fn add_ancestry(
    input_id: IdType,
    left: Position,
    right: Position,
    node: IdType,
    ancestry: &mut AncestryList,
) {
    let head = ancestry.head(input_id).unwrap();
    if head == AncestryList::null() {
        let seg = Segment { left, right, node };
        ancestry.extend(input_id, seg).unwrap();
    } else {
        let last_idx = ancestry.tail(input_id).unwrap();
        if last_idx == AncestryList::null() {
            panic!("last_idx is NULL_ID");
        }
        let last = ancestry.fetch_mut(last_idx).unwrap();
        if last.right == left && last.node == node {
            last.right = right;
        } else {
            let seg = Segment { left, right, node };
            ancestry.extend(input_id, seg).unwrap();
        }
    }
}

fn buffer_edge(
    left: Position,
    right: Position,
    parent: IdType,
    child: IdType,
    temp_edge_buffer: &mut EdgeTable,
) {
    let i = temp_edge_buffer
        .iter()
        .rposition(|e: &Edge| e.child == child);

    match i {
        None => temp_edge_buffer.push(Edge {
            left,
            right,
            parent,
            child,
        }),
        Some(x) => {
            if temp_edge_buffer[x].right == left {
                temp_edge_buffer[x].right = right;
            } else {
                temp_edge_buffer.push(Edge {
                    left,
                    right,
                    parent,
                    child,
                });
            }
        }
    }
}

fn output_buffered_edges(temp_edge_buffer: &mut EdgeTable, new_edges: &mut EdgeTable) -> usize {
    temp_edge_buffer.sort_by(|a, b| a.child.cmp(&b.child));

    // Need to store size here b/c
    // append drains contents of input!!!
    let rv = temp_edge_buffer.len();
    new_edges.append(temp_edge_buffer);

    rv
}

pub fn merge_ancestors(
    input_nodes: &[Node],
    maxlen: Position,
    parent_input_id: IdType,
    state: &mut SimplificationBuffers,
    idmap: &mut [IdType],
) {
    let mut output_id = idmap[parent_input_id as usize];
    let is_sample = output_id != NULL_ID;

    if is_sample {
        state.ancestry.nullify_list(parent_input_id).unwrap();
    }

    let mut previous_right: Position = 0;
    let mut ancestry_node: IdType;
    state.overlapper.init();
    state.temp_edge_buffer.clear();

    while state.overlapper.advance() {
        if state.overlapper.num_overlaps() == 1 {
            ancestry_node = state.overlapper.overlap(0).node;
            if is_sample {
                buffer_edge(
                    state.overlapper.get_left(),
                    state.overlapper.get_right(),
                    output_id,
                    ancestry_node,
                    &mut state.temp_edge_buffer,
                );
                ancestry_node = output_id;
            }
        } else {
            if output_id == NULL_ID {
                state.new_nodes.push(Node {
                    time: input_nodes[parent_input_id as usize].time,
                    deme: input_nodes[parent_input_id as usize].deme,
                });
                output_id = (state.new_nodes.len() - 1) as IdType;
                idmap[parent_input_id as usize] = output_id;
            }
            ancestry_node = output_id;
            for i in 0..state.overlapper.num_overlaps() as usize {
                let o = &state.overlapper.overlap(i);
                buffer_edge(
                    state.overlapper.get_left(),
                    state.overlapper.get_right(),
                    output_id,
                    o.node,
                    &mut state.temp_edge_buffer,
                );
            }
        }
        if is_sample && state.overlapper.get_left() != previous_right {
            add_ancestry(
                parent_input_id,
                previous_right,
                state.overlapper.get_left(),
                output_id,
                &mut state.ancestry,
            );
        }
        add_ancestry(
            parent_input_id,
            state.overlapper.get_left(),
            state.overlapper.get_right(),
            ancestry_node,
            &mut state.ancestry,
        );
        previous_right = state.overlapper.get_right();
    }
    if is_sample && previous_right != maxlen {
        add_ancestry(
            parent_input_id,
            previous_right,
            maxlen,
            output_id,
            &mut state.ancestry,
        );
    }

    if output_id != NULL_ID {
        let n = output_buffered_edges(&mut state.temp_edge_buffer, &mut state.new_edges);

        if n == 0 && !is_sample {
            assert!(output_id < state.new_nodes.len() as IdType);
            state.new_nodes.truncate(output_id as usize);
            idmap[parent_input_id as usize] = NULL_ID;
        }
    }
}

// FIXME: we are panic!-ing here, too
pub fn record_sample_nodes(
    samples: &[IdType],
    tables: &TableCollection,
    new_nodes: &mut NodeTable,
    ancestry: &mut AncestryList,
    idmap: &mut [IdType],
) {
    for sample in samples.iter() {
        assert!(*sample >= 0);
        // NOTE: the following can be debug_assert?
        if *sample == NULL_ID {
            panic!("sample node is NULL_ID");
        }
        if idmap[*sample as usize] != NULL_ID {
            panic!("invalid sample list!");
        }
        let n = tables.node(*sample);
        new_nodes.push(Node {
            time: n.time,
            deme: n.deme,
        });

        add_ancestry(
            *sample,
            0,
            tables.get_length(),
            (new_nodes.len() - 1) as IdType,
            ancestry,
        );

        idmap[*sample as usize] = (new_nodes.len() - 1) as IdType;
    }
}
