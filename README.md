# Forward simulation with tree sequence recording in rust

This package is currently "experimental"!

`forrustts` (pronounced "forests") is a port of many ideas from the [fwdpp](https://github.com/ForwardSimulation/fwdpp) library from C++ to rust.

It is licensed under the MIT license.

## Packaging

* [crate](https://crates.io/crates/forrustts)
* [docs](https://docs.rs/forrustts)

## Getting started

```
cargo build
cargo test --workspace
```

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

