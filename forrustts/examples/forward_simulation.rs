#[path = "../tests/neutral_wright_fisher.rs"]
mod neutral_wright_fisher;

use clap::{Arg, Command};
use neutral_wright_fisher::*;

fn main() {
    let matches = Command::new("forward_simulation")
        .arg(
            Arg::new("popsize")
                .short('N')
                .long("popsize")
                .help("Diploid population size")
                .takes_value(true),
        )
        .arg(
            Arg::new("nsteps")
                .short('n')
                .long("nsteps")
                .help("number of birth steps to simulate")
                .takes_value(true),
        )
        .arg(
            Arg::new("xovers")
                .short('x')
                .long("xovers")
                .help("Mean number of crossovers per meiosis.")
                .takes_value(true),
        )
        .arg(
            Arg::new("mutrate")
                .short('m')
                .long("mutrate")
                .help("Mean number of mutations per new gamete.")
                .takes_value(true),
        )
        .arg(
            Arg::new("simplification_interval")
                .short('s')
                .long("simplify")
                .help("number of generations between simplifications")
                .takes_value(true),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .help("Name of output file. The format is a tskit \"trees\" file")
                .takes_value(true),
        )
        .arg(
            Arg::new("seed")
                .short('S')
                .long("seed")
                .help("Random number seed")
                .takes_value(true),
        )
        .arg(
            Arg::new("psurvival")
                .short('P')
                .long("psurvival")
                .help("Survival probability")
                .takes_value(true),
        )
        .arg(
            Arg::new("validate_tables")
                .short('v')
                .long("validate_tables")
                .help("Validate all tables prior to simplification")
                .takes_value(false),
        )
        .get_matches();

    // TODO: default params

    let popsize = matches.value_of_t("popsize").unwrap();
    let nsteps = matches.value_of_t("nsteps").unwrap();
    let xovers = matches.value_of_t("xovers").unwrap();
    let mutrate = matches.value_of_t("mutrate").unwrap();
    let simplify_input = matches.value_of_t("simplification_interval").unwrap_or(100);
    let psurvival = matches.value_of_t("psurvival").unwrap_or(0.0);
    let seed = matches.value_of_t("seed").unwrap();
    let outfile: String = matches.value_of_t("outfile").unwrap();
    let validate_tables = matches.is_present("validate_tables");

    // TODO: parameter validation..

    let simplify = match simplify_input > 0 {
        true => Some(simplify_input as i64),
        false => None,
    };

    if simplify.is_none() {
        panic!(
            "Simplification interval must be > 0. We got {}",
            simplify_input
        );
    }

    let mut simplification_flags = forrustts_tables_trees::SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= forrustts_tables_trees::SimplificationFlags::VALIDATE_ALL;
    }

    let simparams = SimulationParams {
        popsize,
        mutrate,
        psurvival,
        xovers,
        genome_length: 10000000.into(),
        buffer_edges: true,
        simplification_interval: simplify,
        seed,
        nsteps,
        flags: SimulationFlags::USE_STATE | SimulationFlags::BUFFER_EDGES,
        simplification_flags: forrustts_tables_trees::SimplificationFlags::empty(),
    };

    let (tables, is_sample) = neutral_wf(simparams).unwrap();

    let tskit_tables = forrustts_tskit::export_tables(
        tables,
        forrustts_tskit::simple_time_reverser(nsteps),
        match simplify.is_some() {
            true => Some(forrustts_tskit::TableCollectionExportFlags::BUILD_INDEXES),
            false => None,
        },
    )
    .unwrap();

    let ts = match tskit_tables.tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES) {
        Ok(x) => x,
        Err(e) => panic!("{}", e),
    };

    for i in ts.sample_nodes().iter() {
        assert!(is_sample[usize::from(*i)] > 0);
    }

    ts.dump(&outfile, tskit::TableOutputOptions::default())
        .unwrap();
}
