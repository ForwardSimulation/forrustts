use crate::stochastic_testing_tools::*;
use crate::*;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::Uniform;
use streaming_iterator::StreamingIterator;
use tskit::TableAccess;
use tskit::TskitTypeAccess;

fn compare_edge_table_indexes(
    tables: &TableCollection,
    tsk_tables: &tskit::TableCollection,
) -> bool {
    let tsk_edge_input = unsafe {
        std::slice::from_raw_parts(
            (*tsk_tables.as_ptr()).indexes.edge_insertion_order,
            tsk_tables.edges().num_rows() as usize,
        )
    };
    let tsk_edge_output = unsafe {
        std::slice::from_raw_parts(
            (*tsk_tables.as_ptr()).indexes.edge_removal_order,
            tsk_tables.edges().num_rows() as usize,
        )
    };
    for (idx, val) in tables.edge_input_order.iter().enumerate() {
        assert_eq!(
            tables.edges_[*val].parent,
            tsk_tables.edges().parent(tsk_edge_input[idx]).unwrap()
        );
        assert_eq!(
            tables.edges_[*val].child,
            tsk_tables.edges().child(tsk_edge_input[idx]).unwrap()
        );
        assert_eq!(
            tables.edges_[*val].left,
            tsk_tables.edges().left(tsk_edge_input[idx]).unwrap() as Position
        );
        assert_eq!(
            tables.edges_[*val].right,
            tsk_tables.edges().right(tsk_edge_input[idx]).unwrap() as Position
        );
    }

    for (idx, val) in tables.edge_output_order.iter().enumerate() {
        assert_eq!(
            tables.edges_[*val].parent,
            tsk_tables.edges().parent(tsk_edge_output[idx]).unwrap()
        );
        assert_eq!(
            tables.edges_[*val].child,
            tsk_tables.edges().child(tsk_edge_output[idx]).unwrap()
        );
        assert_eq!(
            tables.edges_[*val].left,
            tsk_tables.edges().left(tsk_edge_output[idx]).unwrap() as Position
        );
        assert_eq!(
            tables.edges_[*val].right,
            tsk_tables.edges().right(tsk_edge_output[idx]).unwrap() as Position
        );
    }
    true
}

#[test]
#[ignore]
fn compare_state_to_no_state() {
    let params = SimulationParams {
        popsize: 1000,
        mutrate: 0.,
        psurvival: 0.5,
        xovers: 5e-3,
        genome_length: 1000000,
        buffer_edges: false,
        simplification_interval: Some(101),
        seed: 666,
        nsteps: 2000,
        flags: SimulationFlags::empty(),
        simplification_flags: SimplificationFlags::empty(),
    };

    let mut params_state = params;
    params_state.flags = SimulationFlags::USE_STATE;

    let (tables, _, _) = neutral_wf(params).unwrap();
    let (tables_state, _, _) = neutral_wf(params_state).unwrap();
    assert_eq!(tables.num_nodes(), tables_state.num_nodes());
    assert_eq!(tables.num_edges(), tables_state.num_edges());

    for (i, j) in tables.nodes().iter().zip(tables_state.nodes()) {
        assert_eq!(i.time, j.time);
        assert_eq!(i.deme, j.deme);
    }
    for (i, j) in tables.edges().iter().zip(tables_state.edges()) {
        assert_eq!(i.left, j.left);
        assert_eq!(i.right, j.right);
        assert_eq!(i.parent, j.parent);
        assert_eq!(i.child, j.child);
    }
}

#[test]
#[ignore]
fn compare_buffer_vs_sort_overlapping_gens() {
    let params = SimulationParams {
        popsize: 1000,
        mutrate: 0.,
        psurvival: 0.5,
        xovers: 5e-3,
        genome_length: 1000000,
        buffer_edges: false,
        simplification_interval: Some(101),
        seed: 666,
        nsteps: 2000,
        flags: SimulationFlags::USE_STATE,
        simplification_flags: SimplificationFlags::empty(),
    };

    let mut params_buffer = params;
    params_buffer.flags = SimulationFlags::BUFFER_EDGES;
    let (mut tables, _, _) = neutral_wf(params).unwrap();
    let (mut tables_buffer, _, _) = neutral_wf(params_buffer).unwrap();
    assert_eq!(tables.num_nodes(), tables_buffer.num_nodes());
    assert_eq!(tables.num_edges(), tables_buffer.num_edges());

    tables.build_indexes(IndexTablesFlags::empty()).unwrap();
    tables_buffer
        .build_indexes(IndexTablesFlags::empty())
        .unwrap();
    assert_eq!(tables.count_trees(), tables_buffer.count_trees());

    let ts = TreeSequence::new(&tables, TreeSequenceFlags::empty()).unwrap();
    let ts_buffer = TreeSequence::new(&tables_buffer, TreeSequenceFlags::empty()).unwrap();

    let mut ti = ts.tree_iterator(TreeFlags::empty());
    let mut ti_buffer = ts_buffer.tree_iterator(TreeFlags::empty());

    for _ in 0..ts.num_trees() {
        if let Some(tree) = ti.next() {
            if let Some(tree_buffer) = ti_buffer.next() {
                assert_eq!(
                    tree.total_branch_length(false).unwrap(),
                    tree_buffer.total_branch_length(false).unwrap()
                );
            } else {
                panic!("expected a Tree");
            }
        } else {
            panic!("expected a Tree");
        }
    }
}

#[test]
#[ignore]
fn simplify_to_samples() {
    let params = SimulationParams {
        popsize: 250,
        mutrate: 2e-3,
        psurvival: 0.5,
        xovers: 5e-3,
        genome_length: 1000000,
        buffer_edges: false,
        simplification_interval: None,
        seed: 1512512,
        nsteps: 500,
        flags: SimulationFlags::empty(),
        simplification_flags: SimplificationFlags::empty(),
    };

    let sims = Simulator::new(params, 10);

    for mut i in sims.iter() {
        let samples = make_samples(&i.is_sample);
        assert!(!samples.samples.is_empty());
        let mut output = SimplificationOutput::new();
        simplify_tables_without_state(
            &samples,
            SimplificationFlags::empty(),
            &mut i.tables,
            &mut output,
        )
        .unwrap();

        i.tsk_tables
            .simplify(
                &samples.samples,
                tskit::SimplificationOptions::FILTER_SITES,
                false,
            )
            .unwrap();

        assert_eq!(
            i.tables.edges().len(),
            i.tsk_tables.edges().num_rows() as usize
        );
        assert_eq!(
            i.tables.nodes().len(),
            i.tsk_tables.nodes().num_rows() as usize
        );
        assert_eq!(
            i.tables.sites().len(),
            i.tsk_tables.sites().num_rows() as usize
        );
        assert_eq!(
            i.tables.mutations().len(),
            i.tsk_tables.mutations().num_rows() as usize
        );
        for (idx, s) in i.tables.enumerate_sites() {
            match i
                .tsk_tables
                .sites()
                .position(idx as tskit::tsk_id_t)
                .unwrap()
                .partial_cmp(&(s.position as f64))
            {
                None => panic!("bad cmp"),
                Some(std::cmp::Ordering::Equal) => (),
                Some(_) => panic!("Expected Equal"),
            };
        }

        for (idx, m) in i.tables.enumerate_mutations() {
            assert_eq!(
                m.node,
                i.tsk_tables
                    .mutations()
                    .node(idx as tskit::tsk_id_t)
                    .unwrap()
            );
        }

        for (idx, m) in i.tables.enumerate_mutations() {
            let tpos = i
                .tsk_tables
                .sites()
                .position(
                    i.tsk_tables
                        .mutations()
                        .site(idx as tskit::tsk_id_t)
                        .unwrap(),
                )
                .unwrap();
            match tpos.partial_cmp(&(i.tables.site(m.site as IdType).position as f64)) {
                Some(std::cmp::Ordering::Equal) => (),
                Some(_) => panic!("Expected Equal"),
                None => panic!("Expected Equal"),
            }
        }

        // Detailed comparisons
        i.tables.build_indexes(IndexTablesFlags::empty()).unwrap();
        i.tsk_tables.build_index().unwrap();
        assert!(compare_edge_table_indexes(&i.tables, &i.tsk_tables));
        let ts = TreeSequence::new(&i.tables, TreeSequenceFlags::empty()).unwrap();
        let tsk_ts =
            tskit::TreeSequence::new(i.tsk_tables, tskit::TreeSequenceFlags::empty()).unwrap();
        assert_eq!(ts.num_trees(), tsk_ts.num_trees());

        let mut ti = ts.tree_iterator(TreeFlags::TRACK_SAMPLES);
        let mut tsk_ti = tsk_ts
            .tree_iterator(tskit::TreeFlags::SAMPLE_LISTS)
            .unwrap();

        for _ in 0..ts.num_trees() {
            let tree = ti.next().unwrap();
            let tsk_tree = tsk_ti.next().unwrap();

            for s in ts.sample_nodes() {
                let mut p = vec![];
                let mut tsk_p = vec![];
                for u in tree.parents(*s).unwrap() {
                    p.push(u);
                }
                for u in tsk_tree.path_to_root(*s).unwrap() {
                    tsk_p.push(u);
                }
                assert!(p == tsk_p);
                for pi in &p {
                    let mut ci = vec![];
                    let mut tsk_ci = vec![];
                    for child in tree.children(*pi).unwrap() {
                        ci.push(child);
                    }
                    for child in tsk_tree.children(*pi).unwrap() {
                        tsk_ci.push(child);
                    }
                    ci.sort_unstable();
                    tsk_ci.sort_unstable();
                    assert!(ci == tsk_ci);
                }
            }
            for node in tree.traverse_nodes(NodeTraversalOrder::Preorder).take(5) {
                let mut samples = vec![];
                let mut tsk_samples = vec![];
                for s in tree.samples(node).unwrap() {
                    samples.push(s);
                }
                for s in tsk_tree.samples(node).unwrap() {
                    tsk_samples.push(s);
                }
                samples.sort_unstable();
                tsk_samples.sort_unstable();
                assert!(samples == tsk_samples);
            }
        }
    }
}

#[test]
#[ignore]
fn simplify_to_arbitrary_nodes() {
    let params = SimulationParams {
        popsize: 250,
        mutrate: 2e-3,
        psurvival: 0.5,
        xovers: 5e-3,
        genome_length: 1000000,
        buffer_edges: false,
        simplification_interval: None,
        seed: 5312851,
        nsteps: 500,
        flags: SimulationFlags::empty(),
        simplification_flags: SimplificationFlags::empty(),
    };

    let sims = Simulator::new(params, 25);

    let mut rng = StdRng::seed_from_u64(588512852);

    let subsample_size = 50;

    for i in sims.iter() {
        for _ in 0..10 {
            let mut tables = i.tables.clone();
            let mut tsk_tables = i.tsk_tables.deepcopy().unwrap();

            let mut candidate_sample: Vec<IdType> = vec![];
            for (idx, val) in i.is_sample.iter().enumerate() {
                if *val == 1 {
                    candidate_sample.push(idx as IdType);
                }
            }
            let node_sampler = Uniform::new(0_usize, candidate_sample.len());
            let mut subsample = vec![0; candidate_sample.len()];
            for _ in 0..subsample_size {
                let mut x = rng.sample(node_sampler);
                while subsample[x] == 1 {
                    x = rng.sample(node_sampler);
                }
                subsample[x] = 1;
            }
            let mut samples_list: Vec<IdType> = vec![];
            for (idx, val) in subsample.iter().enumerate() {
                if *val == 1 {
                    samples_list.push(candidate_sample[idx]);
                }
            }
            assert_eq!(samples_list.len(), subsample_size as usize);
            let mut samples = SamplesInfo::new();
            samples.samples = samples_list;
            let mut output = SimplificationOutput::new();
            simplify_tables_without_state(
                &samples,
                SimplificationFlags::empty(),
                &mut tables,
                &mut output,
            )
            .unwrap();

            tsk_tables
                .simplify(
                    &samples.samples,
                    tskit::SimplificationOptions::FILTER_SITES,
                    false,
                )
                .unwrap();

            assert_eq!(tables.edges().len(), tsk_tables.edges().num_rows() as usize);
            assert_eq!(tables.nodes().len(), tsk_tables.nodes().num_rows() as usize);
            assert_eq!(tables.sites().len(), tsk_tables.sites().num_rows() as usize);
            assert_eq!(
                tables.mutations().len(),
                tsk_tables.mutations().num_rows() as usize
            );
            for (i, s) in tables.enumerate_sites() {
                match tsk_tables
                    .sites()
                    .position(i as tskit::tsk_id_t)
                    .unwrap()
                    .partial_cmp(&(s.position as f64))
                {
                    None => panic!("bad cmp"),
                    Some(std::cmp::Ordering::Equal) => (),
                    Some(_) => panic!("Expected Equal"),
                };
            }

            for (i, m) in tables.enumerate_mutations() {
                assert_eq!(
                    m.node,
                    tsk_tables.mutations().node(i as tskit::tsk_id_t).unwrap()
                );
            }

            for (i, m) in tables.enumerate_mutations() {
                match tsk_tables
                    .sites()
                    .position(tsk_tables.mutations().site(i as tskit::tsk_id_t).unwrap())
                    .unwrap()
                    .partial_cmp(&(tables.site(m.site as IdType).position as f64))
                {
                    Some(std::cmp::Ordering::Equal) => (),
                    Some(_) => panic!("Expected Equal"),
                    None => panic!("Expected Equal"),
                }
            }
        }
    }
}