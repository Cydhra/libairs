# Tree Sequence Inference in Rust

This repository contains an experimental Rust implementation of the core inference algorithm
of [tsinfer](https://github.com/tskit-dev/tsinfer).
The implementation has an experimental parallelization scheme that improves upon tsinfer's performance on high processor
counts (above 32 cores).