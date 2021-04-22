# Forward simulation with tree sequence recording in rust

This package is currently "experimental"!

`forrustts` (pronounced "forests") is a port of many ideas from the [fwdpp](https://github.com/ForwardSimulation/fwdpp) library from C++ to rust.

It is licensed under the GNU General Public License, version 3 or later ("GPL3+").

## Packaging

* [crate](https://crates.io/crates/forrustts)
* [docs](https://docs.rs/forrustts)

## Getting started

```
cargo build
cargo test
```

Example programs are in the subdirectory `examples/`:

```
cargo build --examples
```

The binaries will then be found in `target/debug/examples`.

To build optimized examples:

```
cargo build --release --examples
```

The binaries will then be found in `target/release/examples`.

These programs use [clap](https://crates.io/crates/clap) for command-line options.
Pass ``--help`` to any of them for usage information.

## Development information

## CI

CI testing is done using GitHub actions for both `Linux` and `macOS`.
These actions include using [clippy](https://crates.io/crates/clippy/), which is a very strict code linter.
The actions also check code format using [rustfmt](https://crates.io/crates/rustfmt-nightly).

### Code coverage

Use [tarpaulin](https://docs.rs/crate/cargo-tarpaulin/).
The documentation for that crate is excellent.
The short version is:

```
cargo tarpaulin -o html
```

This command will run the tests and generate a nice `html` report.

## Change log

See [here](https://github.com/ForwardSimulation/forrustts/blob/main/CHANGELOG.md).

