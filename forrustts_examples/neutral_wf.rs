use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::wright_fisher::*;
use forrustts::SimplificationFlags;

fn main() {
    let matches = App::new("neutral_wf")
        .arg(
            Arg::with_name("popsize")
                .short("N")
                .long("popsize")
                .help("Diploid population size")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("nsteps")
                .short("n")
                .long("nsteps")
                .help("number of birth steps to simulate")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("rho")
                .long("rho")
                .help("Scaled rate of crossing over, 4Nr")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("simplification_interval")
                .short("s")
                .long("simplify")
                .help("number of generations between simplifications")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("outfile")
                .short("o")
                .long("outfile")
                .help("Name of output file. The format is a tskit \"trees\" file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("seed")
                .short("S")
                .long("seed")
                .help("Random number seed")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("buffer_edges")
                .short("B")
                .long("buffer_edges")
                .help("Use edge buffering instead of sorting")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("psurvival")
                .short("P")
                .long("psurvival")
                .help("Survival probability")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("validate_tables")
                .short("v")
                .long("validate_tables")
                .help("Validate all tables prior to simplification")
                .takes_value(false),
        )
        .get_matches();

    // TODO: default params

    let popsize = value_t_or_exit!(matches.value_of("popsize"), u32);
    let g = value_t_or_exit!(matches.value_of("nsteps"), i64);
    let rho = value_t_or_exit!(matches.value_of("rho"), f64);
    let simplify_input = value_t!(matches.value_of("simplification_interval"), i64).unwrap_or(-1);
    let psurvival = value_t!(matches.value_of("psurvival"), f64).unwrap_or(0.0);
    let seed = value_t_or_exit!(matches.value_of("seed"), usize);
    let outfile = value_t_or_exit!(matches.value_of("outfile"), String);
    let validate_tables = matches.is_present("validate_tables");

    // TODO: parameter validation..

    let mut simplify: Option<i64> = None;

    if simplify_input > 0 {
        simplify = Some(simplify_input);
    }

    let r = rho / (4.0 * popsize as f64);

    let flags = if matches.is_present("buffer_edges") {
        SimulationFlags::BUFFER_EDGES
    } else {
        SimulationFlags::USE_STATE
    };

    let mut simplification_flags = SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= SimplificationFlags::VALIDATE_ALL;
    }

    let (mut tables, is_sample) = neutral_wf(
        PopulationParams {
            size: popsize,
            genome_length: 10000000,
            littler: r,
            psurvival,
        },
        SimulationParams {
            simplification_interval: simplify,
            seed,
            nsteps: g,
            flags,
            simplification_flags,
        },
    )
    .unwrap();

    let mut tskit_tables = forrustts::tskit_tools::convert_to_tskit_and_drain_minimal(
        &is_sample,
        forrustts::tskit_tools::simple_time_reverser(g),
        simplify.is_some(),
        &mut tables,
    );

    match simplify.is_some() {
        true => (),
        false => {
            let _ = tskit_tables.build_index().unwrap();
        }
    };

    tskit_tables
        .dump(&outfile, tskit::TableOutputOptions::default())
        .unwrap();
}
