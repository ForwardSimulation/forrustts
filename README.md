# Forward simulation with tree sequence recording in rust

This package is currently "experimental"!

`forrustts` (pronounced "forests") is a port of many ideas from the [fwdpp](https://github.com/molpopgen/fwdpp) library from C++ to rust.

It is licensed under the GNU General Public License, version 3 or later ("GPL3+").

## Getting started

Install the GNU Scientific Library.
For example:

```
apt install libgsl-dev
```

You may use `conda` or `brew` as you see fit.

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

## Change log

### 2021-03-26, Version 0.1.2

#### Commits
- [[`21a77982a8`](https://github.com/molpopgen/forrustts/commit/21a77982a8530ae32791ea7eee35af8b9f6da0a7)] bump version to 0.1.2 (molpopogen)
- [[`4653a8b5d1`](https://github.com/molpopgen/forrustts/commit/4653a8b5d128405e9da6f50c58caadb92d97beca)] Merge pull request #64 from molpopgen/update_tskit_dependency (Kevin R. Thornton)
- [[`6bc5030d9e`](https://github.com/molpopgen/forrustts/commit/6bc5030d9e349a0d059f1d06cb741a94522cee01)] Fix string format lint in assert. (molpopogen)
- [[`f36c6c91df`](https://github.com/molpopgen/forrustts/commit/f36c6c91df6500ecc6f2c31246bcf5d5b32bd59e)] update to tskit 0.1.1 (molpopogen)
- [[`b1424be2f8`](https://github.com/molpopgen/forrustts/commit/b1424be2f80400c799503c08500cf48b86c6cfd1)] Merge pull request #62 from molpopgen/doc_fixes (Kevin R. Thornton)
- [[`5771eb18d6`](https://github.com/molpopgen/forrustts/commit/5771eb18d6659af55912697bb7bfeac2b76f2c48)] Document other differences from tskit. (molpopgen)


