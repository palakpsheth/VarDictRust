#![no_main]

use std::io::Cursor;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::bam::{RecordPreprocessor, SamRecord};
use vardictrust_core::config::Config;

fuzz_target!(|data: &[u8]| {
    let header_src = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n";
    let mut sam_bytes = Vec::with_capacity(header_src.len() + data.len() + 1);
    sam_bytes.extend_from_slice(header_src);
    sam_bytes.extend_from_slice(data);
    sam_bytes.push(b'\n');

    let mut reader = noodles::sam::io::Reader::new(Cursor::new(sam_bytes));
    let Ok(header) = reader.read_header() else {
        return;
    };

    let config = Config::default();
    let mut preprocessor = RecordPreprocessor::new();
    let mut sam_rec = SamRecord::new();

    for result in reader.records().take(16) {
        if let Ok(record) = result {
            if sam_rec.fill_from_alignment_record(&record, &header).is_ok() {
                let _ = preprocessor.preprocess_record(&sam_rec, &config);
            }
        }
    }
});
