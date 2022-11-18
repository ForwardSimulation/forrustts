#[path = "./stochastic_testing_tools.rs"]
mod stochastic_testing_tools;

use std::convert::TryFrom;

use forrustts::NodeId;
use forrustts::*;
use rand::Rng;
use rand::SeedableRng;
use stochastic_testing_tools::*;
use streaming_iterator::StreamingIterator;
use tskit::bindings::tsk_id_t;
use tskit::TskitTypeAccess;

struct VecTskId(Vec<tsk_id_t>);

impl From<Vec<NodeId>> for VecTskId {
    fn from(value: Vec<NodeId>) -> Self {
        let mut rv: Vec<tsk_id_t> = vec![];
        for v in &value {
            rv.push(v.raw());
        }
        Self(rv)
    }
}

fn compare_edge_table_indexes(
    tables: &TableCollection,
    tsk_tables: &tskit::TableCollection,
) -> bool {
    let tsk_edge_input = unsafe {
        std::slice::from_raw_parts(
            (*tsk_tables.as_ptr()).indexes.edge_insertion_order,
            usize::try_from(tsk_tables.edges().num_rows()).unwrap(),
        )
    };
    let tsk_edge_output = unsafe {
        std::slice::from_raw_parts(
            (*tsk_tables.as_ptr()).indexes.edge_removal_order,
            usize::try_from(tsk_tables.edges().num_rows()).unwrap(),
        )
    };
    for (idx, val) in tables.edge_input_order().unwrap().iter().enumerate() {
        assert_eq!(
            tables.edges()[*val].parent.raw(),
            tsk_tables.edges().parent(tsk_edge_input[idx]).unwrap()
        );
        assert_eq!(
            tables.edges()[*val].child.raw(),
            tsk_tables.edges().child(tsk_edge_input[idx]).unwrap()
        );
        assert_eq!(
            tables.edges()[*val].left.raw(),
            f64::from(tsk_tables.edges().left(tsk_edge_input[idx]).unwrap()) as i64
        );
        assert_eq!(
            tables.edges()[*val].right.raw(),
            f64::from(tsk_tables.edges().right(tsk_edge_input[idx]).unwrap()) as i64
        );
    }

    for (idx, val) in tables.edge_output_order().unwrap().iter().enumerate() {
        assert_eq!(
            tables.edges()[*val].parent.raw(),
            tsk_tables.edges().parent(tsk_edge_output[idx]).unwrap()
        );
        assert_eq!(
            tables.edges()[*val].child.raw(),
            tsk_tables.edges().child(tsk_edge_output[idx]).unwrap()
        );
        assert_eq!(
            tables.edges()[*val].left.raw(),
            f64::from(tsk_tables.edges().left(tsk_edge_output[idx]).unwrap()) as i64
        );
        assert_eq!(
            tables.edges()[*val].right.raw(),
            f64::from(tsk_tables.edges().right(tsk_edge_output[idx]).unwrap()) as i64
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
        genome_length: 1000000.into(),
        buffer_edges: false,
        simplification_interval: Some(101),
        seed: 666,
        nsteps: 2000,
        flags: SimulationFlags::empty(),
        simplification_flags: SimplificationFlags::empty(),
    };

    let mut params_state = params;
    params_state.flags = SimulationFlags::USE_STATE;

    let (tables, _) = neutral_wf(params).unwrap();
    let (tables_state, _) = neutral_wf(params_state).unwrap();
    assert_eq!(tables.num_nodes(), tables_state.num_nodes());
    assert_eq!(tables.num_edges(), tables_state.num_edges());

    for (i, j) in tables.nodes().iter().zip(tables_state.nodes()) {
        match i.time.partial_cmp(&j.time) {
            Some(std::cmp::Ordering::Equal) => (),
            Some(_) => panic!("expected Equal, got not equal"),
            None => panic!("expected Equal, got not None"),
        };
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
        genome_length: 1000000.into(),
        buffer_edges: false,
        simplification_interval: Some(101),
        seed: 666,
        nsteps: 2000,
        flags: SimulationFlags::USE_STATE,
        simplification_flags: SimplificationFlags::empty(),
    };

    let mut params_buffer = params;
    params_buffer.flags = SimulationFlags::BUFFER_EDGES;
    let (mut tables, _) = neutral_wf(params).unwrap();
    let (mut tables_buffer, _) = neutral_wf(params_buffer).unwrap();
    assert_eq!(tables.num_nodes(), tables_buffer.num_nodes());
    assert_eq!(tables.num_edges(), tables_buffer.num_edges());

    tables.build_indexes(IndexTablesFlags::empty()).unwrap();
    tables_buffer
        .build_indexes(IndexTablesFlags::empty())
        .unwrap();
    assert_eq!(tables.count_trees(), tables_buffer.count_trees());

    let ts = TreeSequence::new(tables, TreeSequenceFlags::empty()).unwrap();
    let ts_buffer = TreeSequence::new(tables_buffer, TreeSequenceFlags::empty()).unwrap();

    let mut ti = ts.tree_iterator(TreeFlags::empty());
    let mut ti_buffer = ts_buffer.tree_iterator(TreeFlags::empty());

    for _ in 0..ts.num_trees() {
        if let Some(tree) = ti.next() {
            if let Some(tree_buffer) = ti_buffer.next() {
                let time1 = tree.total_branch_length(false).unwrap();
                let time2 = tree_buffer.total_branch_length(false).unwrap();
                assert!((f64::from(time1) - f64::from(time2)).abs() < f64::EPSILON);
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
        genome_length: 1000000.into(),
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

        let num_input_mutations = i.tables.mutations().len();
        simplify_tables_without_state(
            &samples,
            SimplificationFlags::empty(),
            &mut i.tables,
            &mut output,
        )
        .unwrap();

        assert_eq!(
            output.extinct_mutations.len(),
            num_input_mutations - i.tables.mutations().len()
        );

        i.tsk_tables
            .simplify(
                &VecTskId::from(samples.samples).0,
                tskit::SimplificationOptions::FILTER_SITES,
                false,
            )
            .unwrap();

        assert_eq!(
            i.tables.edges().len(),
            usize::try_from(i.tsk_tables.edges().num_rows()).unwrap()
        );
        assert_eq!(
            i.tables.nodes().len(),
            usize::try_from(i.tsk_tables.nodes().num_rows()).unwrap()
        );
        assert_eq!(
            i.tables.sites().len(),
            usize::try_from(i.tsk_tables.sites().num_rows()).unwrap()
        );
        assert_eq!(
            i.tables.mutations().len(),
            usize::try_from(i.tsk_tables.mutations().num_rows()).unwrap()
        );
        for (idx, s) in i.tables.enumerate_sites() {
            match i
                .tsk_tables
                .sites()
                .position(idx as tsk_id_t)
                .unwrap()
                .partial_cmp(&(s.position.raw() as f64))
            {
                None => panic!("bad cmp"),
                Some(std::cmp::Ordering::Equal) => (),
                Some(_) => panic!("Expected Equal"),
            };
        }

        for (idx, m) in i.tables.enumerate_mutations() {
            assert_eq!(
                m.node.raw(),
                i.tsk_tables.mutations().node(idx as tsk_id_t).unwrap()
            );
        }

        for (idx, m) in i.tables.enumerate_mutations() {
            let tpos = i
                .tsk_tables
                .sites()
                .position(i.tsk_tables.mutations().site(idx as tsk_id_t).unwrap())
                .unwrap();
            match tpos.partial_cmp(&(i.tables.site(m.site).position.raw() as f64)) {
                Some(std::cmp::Ordering::Equal) => (),
                Some(_) => panic!("Expected Equal"),
                None => panic!("Expected Equal"),
            }
        }

        // Detailed comparisons
        i.tables.build_indexes(IndexTablesFlags::empty()).unwrap();
        i.tsk_tables.build_index().unwrap();
        assert!(compare_edge_table_indexes(&i.tables, &i.tsk_tables));
        let ts = TreeSequence::new(i.tables, TreeSequenceFlags::empty()).unwrap();
        let tsk_ts =
            tskit::TreeSequence::new(i.tsk_tables, tskit::TreeSequenceFlags::empty()).unwrap();
        assert_eq!(
            ts.num_trees(),
            usize::try_from(tsk_ts.num_trees()).unwrap() as u32
        );

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
                for u in tree.parents(*s) {
                    p.push(u);
                }
                for u in tsk_tree.parents(s.raw().into()).unwrap() {
                    tsk_p.push(u);
                }
                for (i, j) in p.iter().zip(tsk_p.iter()) {
                    assert_eq!(i.raw(), *j);
                }
                for pi in &p {
                    let mut ci = vec![];
                    let mut tsk_ci = vec![];
                    for child in tree.children(*pi) {
                        ci.push(child);
                    }
                    for child in tsk_tree.children(pi.raw().into()).unwrap() {
                        tsk_ci.push(child);
                    }
                    ci.sort_unstable();
                    tsk_ci.sort_unstable();
                    for (i, j) in ci.iter().zip(tsk_ci.iter()) {
                        assert!(i.raw() == *j);
                    }
                }
            }
            for node in tree.traverse_nodes(NodeTraversalOrder::Preorder).take(5) {
                let mut samples = vec![];
                let mut tsk_samples = vec![];
                for s in tree.samples(node).unwrap() {
                    samples.push(s);
                }
                if let Some(Ok(s)) = tsk_tree.samples(node.raw().into()) {
                    for i in s {
                        tsk_samples.push(i);
                    }
                }
                samples.sort_unstable();
                tsk_samples.sort_unstable();
                for (i, j) in samples.iter().zip(tsk_samples.iter()) {
                    assert_eq!(i.raw(), *j);
                }
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
        genome_length: 1000000.into(),
        buffer_edges: false,
        simplification_interval: None,
        seed: 5312851,
        nsteps: 500,
        flags: SimulationFlags::empty(),
        simplification_flags: SimplificationFlags::empty(),
    };

    let sims = Simulator::new(params, 25);

    let mut rng = rand::rngs::StdRng::seed_from_u64(588512852);

    let subsample_size = 50;

    for i in sims.iter() {
        for _ in 0..10 {
            let mut tables = i.tables.clone();
            let mut tsk_tables = i.tsk_tables.deepcopy().unwrap();

            let mut candidate_sample: Vec<NodeId> = vec![];
            for (idx, val) in i.is_sample.iter().enumerate() {
                if *val == 1 {
                    candidate_sample.push(NodeId::try_from(idx).unwrap());
                }
            }
            //let node_sampler = Uniform::new(0_usize, candidate_sample.len());
            let mut subsample = vec![0; candidate_sample.len()];
            let flat = rand::distributions::Uniform::new(0, candidate_sample.len());
            for _ in 0..subsample_size {
                //let mut x = rng.sample(node_sampler);
                let mut x = rng.sample(flat);
                while subsample[x] == 1 {
                    x = rng.sample(flat);
                }
                subsample[x] = 1;
            }
            let mut samples_list: Vec<NodeId> = vec![];
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
                    &VecTskId::from(samples.samples).0,
                    tskit::SimplificationOptions::FILTER_SITES,
                    false,
                )
                .unwrap();

            assert_eq!(
                tables.edges().len(),
                usize::try_from(tsk_tables.edges().num_rows()).unwrap()
            );
            assert_eq!(
                tables.nodes().len(),
                usize::try_from(tsk_tables.nodes().num_rows()).unwrap()
            );
            assert_eq!(
                tables.sites().len(),
                usize::try_from(tsk_tables.sites().num_rows()).unwrap()
            );
            assert_eq!(
                tables.mutations().len(),
                usize::try_from(tsk_tables.mutations().num_rows()).unwrap()
            );
            for (i, s) in tables.enumerate_sites() {
                match tsk_tables
                    .sites()
                    .position(i as tsk_id_t)
                    .unwrap()
                    .partial_cmp(&(i64::from(s.position) as f64))
                {
                    None => panic!("bad cmp"),
                    Some(std::cmp::Ordering::Equal) => (),
                    Some(_) => panic!("Expected Equal"),
                };
            }

            for (i, m) in tables.enumerate_mutations() {
                assert_eq!(
                    m.node.raw(),
                    tsk_tables.mutations().node(i as tsk_id_t).unwrap()
                );
            }

            for (i, m) in tables.enumerate_mutations() {
                match tsk_tables
                    .sites()
                    .position(tsk_tables.mutations().site(i as tsk_id_t).unwrap())
                    .unwrap()
                    .partial_cmp(&(tables.site(m.site).position.raw() as f64))
                {
                    Some(std::cmp::Ordering::Equal) => (),
                    Some(_) => panic!("Expected Equal"),
                    None => panic!("Expected Equal"),
                }
            }
        }
    }
}

#[test]
#[ignore]
fn test_mutation_tables() {
    let seeds: Vec<u64> = vec![18822, 6699, 173, 14199, 5046, 32637, 25950];
    for seed in seeds {
        let nsteps = 1000;
        let simparams = SimulationParams {
            popsize: 100,
            mutrate: 2e-3,
            psurvival: 0.0,
            xovers: 3e-2,
            genome_length: 10000000.into(),
            buffer_edges: true,
            simplification_interval: Some(100),
            seed,
            nsteps,
            flags: SimulationFlags::USE_STATE | SimulationFlags::BUFFER_EDGES,
            simplification_flags: forrustts::SimplificationFlags::empty(),
        };
        let (tables, is_sample) = neutral_wf(simparams).unwrap();
        forrustts::validate_site_table(tables.genome_length(), tables.sites()).unwrap_or_else(
            |e| {
                panic!("{}", e);
            },
        );
        forrustts::validate_mutation_table(tables.mutations(), tables.sites(), tables.nodes())
            .unwrap_or_else(|e| {
                panic!("{}", e);
            });
        let tskit_tables = forrustts_tskit::export_tables(
            tables.clone(),
            forrustts_tskit::simple_time_reverser(nsteps),
            forrustts_tskit::TableCollectionExportFlags::BUILD_INDEXES,
        )
        .unwrap();
        assert!(tskit_tables.sites().num_rows() > 0);
        assert!(tskit_tables.mutations().num_rows() > 0);

        let ts = tskit_tables
            .tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES)
            .unwrap_or_else(|e| panic!("{}", e));
        for (i, j) in is_sample.iter().enumerate() {
            let flags = ts.nodes().flags(tskit::NodeId::from(i as i32)).unwrap();
            if *j == 1 {
                assert!(flags.contains(tskit::NodeFlags::IS_SAMPLE));
            } else {
                assert!(!flags.contains(tskit::NodeFlags::IS_SAMPLE));
            }
        }
    }
}
