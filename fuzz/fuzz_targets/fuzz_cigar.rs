//! Fuzz target: `parse_cigar` must return Some/None, never panic.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(s) = std::str::from_utf8(data) {
        let parsed = vardict_core::cigar::parse_cigar(s);
        if let Some(ops) = parsed {
            // Invariant: every op length parses as u32 and the element count
            // can never exceed the input length (one element per op char).
            debug_assert!(ops.len() <= s.len());
        }
    }
});
