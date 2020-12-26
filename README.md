# Forward simulation with tree sequence recording in rust

This package is currently "experimental"!

`forrustts` (pronounced "forests") is a port of many ideas from the [fwdpp](https://github.com/molpopgen/fwdpp) library from C++ to rust.

It is licensed under the GNU General Public License, version 3 or later ("GPL3+").

## Getting started

[Install rust](https://www.rust-lang.org/learn/get-started)

Then:

```
cargo build
cargo test
```

Example programs are in the subdirectory `forrustts_examples`:

```
cd forrustts_examples
cargo build --release
```

The binaries will then be found in `target/release/`.

These programs use [clap](https://crates.io/crates/clap) for command-line options.
Pass ``--help`` to any of them for usage information.

## Development information

## CI

CI testing is done using GitHub actions for both `Linux` and `macOS`.
These actions include using [clippy](https://crates.io/crates/clippy/0.0.211), which is a very strict code linter.
The actions also check code format using [rustfmt](https://crates.io/crates/rustfmt-nightly).

### Code coverage

Use [tarpaulin](https://docs.rs/crate/cargo-tarpaulin/0.3.12).
The documentation for that crate is excellent.
The short version is:

```
cargo tarpaulin -o html
```

This command will run the tests and generate a nice `html` report.


