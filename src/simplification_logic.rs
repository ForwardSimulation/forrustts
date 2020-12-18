use crate::nested_forward_list::NestedForwardList;
use crate::segment::Segment;
use crate::tables::*;

//FIXME:
type TsInt = i32;
type SamplesVec = Vec<TsInt>;
const NULLTSINT: TsInt = -1;

pub struct SegmentOverlapper {
    segment_queue: Vec<Segment>,
    overlapping: Vec<Segment>,
    left: i64,
    right: i64,
    qbeg: usize,
    qend: usize,
    obeg: usize,
    oend: usize,
}

impl SegmentOverlapper {
    fn set_partition(&mut self) -> i64 {
        let mut tright = std::i64::MAX;
        let mut b: usize = 0;

        for i in 0..self.oend {
            if self.overlapping[i].right > self.left {
                self.overlapping[b] = self.overlapping[i];
                tright = std::cmp::min(tright, self.overlapping[b].right);
                b += 1;
            }
        }

        self.oend = b;

        return tright;
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
        return self.oend - self.obeg;
    }

    // Public interface below

    pub const fn new() -> SegmentOverlapper {
        return SegmentOverlapper {
            segment_queue: vec![],
            overlapping: vec![],
            left: 0,
            right: std::i64::MAX,
            qbeg: std::usize::MAX,
            qend: std::usize::MAX,
            obeg: std::usize::MAX,
            oend: std::usize::MAX,
        };
    }

    pub fn init(&mut self) -> () {
        self.qbeg = 0;
        self.qend = self.segment_queue.len() - 1;
        assert!(self.qend < self.segment_queue.len());
        self.obeg = 0;
        self.oend = 0;
        self.overlapping.clear();
    }

    pub fn enqueue(&mut self, left: i64, right: i64, node: TsInt) -> () {
        self.segment_queue.push(Segment {
            left: left,
            right: right,
            node: node,
        });
    }

    pub fn finalize_queue(&mut self, maxlen: i64) -> () {
        self.segment_queue.sort_by(|a, b| {
            return a.left.cmp(&b.left);
        });
        self.segment_queue.push(Segment {
            left: maxlen,
            right: maxlen + 1,
            node: NULLTSINT,
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
            self.right = std::i64::MAX;
            let tright = self.set_partition();
            if self.num_overlaps() > 0 {
                self.right = tright;
                rv = true
            }
        }

        return rv;
    }

    pub fn get_left(&self) -> i64 {
        return self.left;
    }

    pub fn get_right(&self) -> i64 {
        return self.right;
    }

    pub fn clear_queue(&mut self) -> () {
        self.segment_queue.clear();
    }

    pub fn overlap(&self, i: usize) -> &Segment {
        return &self.overlapping[i];
    }
}

pub type AncestryList = NestedForwardList<Segment>;

// FIXME: another panic! room
pub fn find_parent_child_segment_overlap(
    edges: &EdgeTable<TsInt>,
    edge_index: usize,
    num_edges: usize,
    maxlen: i64,
    u: TsInt,
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
                return true;
            })
            .unwrap();
        i += 1;
    }
    overlapper.finalize_queue(maxlen);
    return i;
}

pub fn setup_idmap(nodes: &NodeTable<TsInt>) -> SamplesVec {
    return vec![NULLTSINT; nodes.len()];
}

// FIXME: this will panic! if we get error from AncestryList
fn add_ancestry(
    input_id: TsInt,
    left: i64,
    right: i64,
    node: TsInt,
    ancestry: &mut AncestryList,
) -> () {
    let head = ancestry.head(input_id).unwrap();
    if head == AncestryList::null() {
        let seg = Segment {
            left: left,
            right: right,
            node: node,
        };
        ancestry.extend(input_id, seg).unwrap();
    } else {
        let last_idx = ancestry.tail(input_id).unwrap();
        if last_idx == AncestryList::null() {
            panic!("last_idx is NULLTSINT");
        }
        let last = ancestry.fetch_mut(last_idx).unwrap();
        if last.right == left && last.node == node {
            last.right = right;
        } else {
            let seg = Segment {
                left: left,
                right: right,
                node: node,
            };
            ancestry.extend(input_id, seg).unwrap();
        }
    }
}

fn buffer_edge(
    left: i64,
    right: i64,
    parent: TsInt,
    child: TsInt,
    temp_edge_buffer: &mut EdgeTable<TsInt>,
) {
    let i = temp_edge_buffer
        .iter()
        .rposition(|e: &Edge<i32>| e.child == child);

    match i {
        None => temp_edge_buffer.push(Edge::<i32> {
            left: left,
            right: right,
            parent: parent,
            child: child,
        }),
        Some(x) => {
            if temp_edge_buffer[x].right == left {
                temp_edge_buffer[x].right = right;
            } else {
                temp_edge_buffer.push(Edge::<i32> {
                    left: left,
                    right: right,
                    parent: parent,
                    child: child,
                });
            }
        }
    }
}

fn output_buffered_edges(
    temp_edge_buffer: &mut EdgeTable<TsInt>,
    new_edges: &mut EdgeTable<TsInt>,
) -> usize {
    temp_edge_buffer.sort_by(|a, b| {
        return a.child.cmp(&b.child);
    });

    // Need to store size here b/c
    // append drains contents of input!!!
    let rv = temp_edge_buffer.len();
    new_edges.append(temp_edge_buffer);

    return rv;
}

pub fn merge_ancestors(
    input_nodes: &NodeTable<TsInt>,
    maxlen: i64,
    parent_input_id: TsInt,
    temp_edge_buffer: &mut EdgeTable<TsInt>,
    new_nodes: &mut NodeTable<TsInt>,
    new_edges: &mut EdgeTable<TsInt>,
    ancestry: &mut AncestryList,
    overlapper: &mut SegmentOverlapper,
    idmap: &mut SamplesVec,
) {
    let mut output_id = idmap[parent_input_id as usize];
    let is_sample = output_id != NULLTSINT;

    if is_sample == true {
        ancestry.nullify_list(parent_input_id).unwrap();
    }

    let mut previous_right: i64 = 0;
    let mut ancestry_node: TsInt;
    overlapper.init();
    temp_edge_buffer.clear();

    while overlapper.advance() == true {
        if overlapper.num_overlaps() == 1 {
            ancestry_node = overlapper.overlap(0).node;
            if is_sample == true {
                buffer_edge(
                    overlapper.get_left(),
                    overlapper.get_right(),
                    output_id,
                    ancestry_node,
                    temp_edge_buffer,
                );
                ancestry_node = output_id;
            }
        } else {
            if output_id == NULLTSINT {
                new_nodes.push(Node::<i32> {
                    time: input_nodes[parent_input_id as usize].time,
                    deme: input_nodes[parent_input_id as usize].deme,
                });
                output_id = (new_nodes.len() - 1) as TsInt;
                idmap[parent_input_id as usize] = output_id;
            }
            ancestry_node = output_id;
            for i in 0..overlapper.num_overlaps() as usize {
                let o = &overlapper.overlap(i);
                buffer_edge(
                    overlapper.get_left(),
                    overlapper.get_right(),
                    output_id,
                    o.node,
                    temp_edge_buffer,
                );
            }
        }
        if is_sample == true && overlapper.get_left() != previous_right {
            add_ancestry(
                parent_input_id,
                previous_right,
                overlapper.get_left(),
                output_id,
                ancestry,
            );
        }
        add_ancestry(
            parent_input_id,
            overlapper.get_left(),
            overlapper.get_right(),
            ancestry_node,
            ancestry,
        );
        previous_right = overlapper.get_right();
    }
    if is_sample == true && previous_right != maxlen {
        add_ancestry(parent_input_id, previous_right, maxlen, output_id, ancestry);
    }

    if output_id != NULLTSINT {
        let n = output_buffered_edges(temp_edge_buffer, new_edges);

        if n == 0 && is_sample == false {
            assert!(output_id < new_nodes.len() as TsInt);
            new_nodes.truncate(output_id as usize);
            idmap[parent_input_id as usize] = NULLTSINT;
        }
    }
}

// FIXME: we are panic!-ing here, too
pub fn record_sample_nodes(
    samples: &SamplesVec,
    tables: &TableCollection,
    new_nodes: &mut NodeTable<TsInt>,
    ancestry: &mut AncestryList,
    idmap: &mut SamplesVec,
) -> () {
    for sample in samples.iter() {
        assert!(*sample >= 0);
        // NOTE: the following can be debug_assert?
        if *sample == NULLTSINT {
            panic!("sample node is NULLTSINT");
        }
        if idmap[*sample as usize] != NULLTSINT {
            panic!("invalid sample list!");
        }
        let n = tables.node(*sample as usize);
        new_nodes.push(Node::<i32> {
            time: n.time,
            deme: n.deme,
        });

        add_ancestry(
            *sample,
            0,
            tables.get_length(),
            (new_nodes.len() - 1) as TsInt,
            ancestry,
        );

        idmap[*sample as usize] = (new_nodes.len() - 1) as TsInt;
    }
}
