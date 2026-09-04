//! Library target so other crates can reuse the EDF reader/writer
//! (e.g. the Vitaport VPD batch converter). The `convert_edf` binary
//! is unaffected; it compiles these modules directly via `mod`.

pub mod edf;
