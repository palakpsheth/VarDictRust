use crate::types::{Pos0, Pos1};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::{QualityScores, Sequence};
use noodles::sam::alignment::Record;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamRecord {
    pub name: Vec<u8>,
    pub flags: u16,
    pub reference_name: Option<String>,
    pub alignment_start: Option<Pos1>,
    pub alignment_end: Option<Pos1>,
    pub mapping_quality: Option<u8>,
    pub cigar: Vec<(u32, CigarKind)>,
    pub sequence: Vec<u8>,
    pub quality_scores: Vec<u8>,
    pub mate_reference_name: Option<String>,
    pub mate_alignment_start: Option<Pos1>,
    pub template_length: i32,
    pub nm: Option<i32>,
    pub sa: Option<String>,
    pub mc: Option<String>,
    pub mq: Option<u8>,
}

impl SamRecord {
    pub fn from_alignment_record<R>(
        record: &R,
        header: &noodles::sam::Header,
    ) -> Result<Self, std::io::Error>
    where
        R: Record,
    {
        let name = if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            bytes.to_vec()
        } else {
            Vec::new()
        };

        let flags = record.flags().map(|f| f.bits()).unwrap_or(0);

        let reference_name = record
            .reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            });

        let alignment_start = if let Some(pos) = record.alignment_start() {
            pos.ok().map(|p| Pos1(p.get() as u32))
        } else {
            None
        };

        let alignment_end = if let Some(pos) = record.alignment_end() {
            pos.ok().map(|p| Pos1(p.get() as u32))
        } else {
            None
        };

        let mapping_quality = if (flags & 0x4) != 0 {
            None
        } else {
            Some(
                record
                    .mapping_quality()
                    .and_then(|m| m.ok())
                    .map(|m| m.get())
                    .unwrap_or(255),
            )
        };

        let mut cigar = Vec::new();
        for op in record.cigar().iter().flatten() {
            cigar.push((op.len() as u32, op.kind()));
        }

        let mut sequence = Vec::new();
        for b in record.sequence().iter() {
            sequence.push(b);
        }

        let mut quality_scores = Vec::new();
        for score in record.quality_scores().iter().flatten() {
            quality_scores.push(score + 33);
        }

        let mate_reference_name = record
            .mate_reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            });

        let mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(|p| Pos1(p.get() as u32))
        } else {
            None
        };

        let template_length = record.template_length().unwrap_or(0);

        let mut nm = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'n', b'M',
                ))
        {
            nm = val.as_int().map(|x| x as i32);
        } else if let Some(Ok(val)) = record
            .data()
            .get(&noodles::sam::alignment::record::data::field::Tag::EDIT_DISTANCE)
        {
            nm = val.as_int().map(|x| x as i32);
        }

        let mut sa = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'S', b'A',
                ))
        {
            match val {
                noodles::sam::alignment::record::data::field::Value::String(s) => {
                    sa = Some(String::from_utf8_lossy(s.as_ref()).into_owned());
                }
                _ => {
                    sa = Some(format!("{:?}", val));
                }
            }
        }

        let mut mc = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'C',
                ))
        {
            match val {
                noodles::sam::alignment::record::data::field::Value::String(s) => {
                    mc = Some(String::from_utf8_lossy(s.as_ref()).into_owned());
                }
                _ => {
                    mc = Some(format!("{:?}", val));
                }
            }
        }

        let mut mq = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'Q',
                ))
        {
            mq = val.as_int().map(|x| x as u8);
        }

        Ok(SamRecord {
            name,
            flags,
            reference_name,
            alignment_start,
            alignment_end,
            mapping_quality,
            cigar,
            sequence,
            quality_scores,
            mate_reference_name,
            mate_alignment_start,
            template_length,
            nm,
            sa,
            mc,
            mq,
        })
    }
}

pub enum IndexedSamReader {
    Bam(noodles::bam::io::IndexedReader<noodles::bgzf::Reader<std::fs::File>>),
    Cram(noodles::cram::io::IndexedReader<std::fs::File>),
}

impl IndexedSamReader {
    pub fn query<'a>(
        &'a mut self,
        header: &'a noodles::sam::Header,
        region: &'a noodles::core::Region,
    ) -> Result<Box<dyn Iterator<Item = Result<SamRecord, std::io::Error>> + 'a>, std::io::Error>
    {
        match self {
            IndexedSamReader::Bam(reader) => {
                let query = reader.query(header, region)?;
                let iter = query.map(|record| {
                    let rec = record?;
                    SamRecord::from_alignment_record(&rec, header)
                });
                Ok(Box::new(iter))
            }
            IndexedSamReader::Cram(reader) => {
                let query = reader.query(header, region)?;
                let iter = query.map(|record| {
                    let rec = record?;
                    SamRecord::from_alignment_record(&rec, header)
                });
                Ok(Box::new(iter))
            }
        }
    }
}

pub struct CachedReader {
    pub reader: IndexedSamReader,
    pub header: noodles::sam::Header,
}

thread_local! {
    static SAM_READER_CACHE: RefCell<HashMap<PathBuf, CachedReader>> = RefCell::new(HashMap::new());
}

pub fn open_reader(
    path: &Path,
    config: &crate::config::Config,
) -> Result<CachedReader, crate::types::VardictError> {
    let path_str = path.to_string_lossy();
    if path_str.ends_with(".cram") || path_str.ends_with(".CRAM") {
        let fasta_path = std::path::Path::new(&config.reference);
        let mut fai_path = fasta_path.to_path_buf();
        let mut ext = fai_path
            .extension()
            .map(|e| e.to_string_lossy().into_owned())
            .unwrap_or_default();
        ext.push_str(".fai");
        fai_path.set_extension(ext);

        let fasta_reader = noodles::fasta::IndexedReader::new(
            std::io::BufReader::new(
                std::fs::File::open(fasta_path).map_err(crate::types::VardictError::Io)?,
            ),
            noodles::fasta::fai::read(&fai_path).map_err(crate::types::VardictError::Io)?,
        );
        use noodles::fasta::repository::adapters::IndexedReader as FastaAdapter;
        use noodles::fasta::Repository;
        let repository = Repository::new(FastaAdapter::new(fasta_reader));

        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;
        Ok(CachedReader {
            reader: IndexedSamReader::Cram(reader),
            header,
        })
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;
        Ok(CachedReader {
            reader: IndexedSamReader::Bam(reader),
            header,
        })
    }
}

pub fn query_region(
    bam_path: &Path,
    chr: &str,
    start: Pos0,
    end: Pos0,
    config: &crate::config::Config,
) -> Result<Vec<SamRecord>, crate::types::VardictError> {
    SAM_READER_CACHE.with(
        |cache| -> Result<Vec<SamRecord>, crate::types::VardictError> {
            let mut cache = cache.borrow_mut();
            if !cache.contains_key(bam_path) {
                let reader = open_reader(bam_path, config)?;
                cache.insert(bam_path.to_path_buf(), reader);
            }
            let cached = match cache.get_mut(bam_path) {
                Some(c) => c,
                None => {
                    return Err(crate::types::VardictError::Other(format!(
                        "BAM cache miss for {}",
                        bam_path.display()
                    )))
                }
            };

            let start_pos = start.0 + 1;
            let end_pos = end.0;
            if start_pos > end_pos {
                return Ok(Vec::new());
            }

            let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
            let region: noodles::core::Region = region_str.parse().map_err(|e| {
                crate::types::VardictError::Parse(format!(
                    "failed to parse region '{}': {}",
                    region_str, e
                ))
            })?;

            let iter = cached
                .reader
                .query(&cached.header, &region)
                .map_err(crate::types::VardictError::Io)?;
            let mut records = Vec::new();
            for record in iter {
                records.push(record.map_err(crate::types::VardictError::Io)?);
            }
            Ok(records)
        },
    )
}

// Helper to reproduce getMateReferenceName(record) from VarDict
pub fn get_mate_reference_name(record: &SamRecord) -> String {
    match (&record.reference_name, &record.mate_reference_name) {
        (None, _) | (_, None) => "*".to_string(),
        (Some(ref_name), Some(mate_ref)) => {
            if ref_name == mate_ref {
                "=".to_string()
            } else {
                mate_ref.clone()
            }
        }
    }
}

// Helper to format CIGAR to a string matching getCigarString()
pub fn get_cigar_string(cigar: &[(u32, CigarKind)]) -> String {
    let mut s = String::new();
    for (len, kind) in cigar {
        s.push_str(&len.to_string());
        let c = match kind {
            CigarKind::Match => 'M',
            CigarKind::Insertion => 'I',
            CigarKind::Deletion => 'D',
            CigarKind::Skip => 'N',
            CigarKind::SoftClip => 'S',
            CigarKind::HardClip => 'H',
            CigarKind::Pad => 'P',
            CigarKind::SequenceMatch => '=',
            CigarKind::SequenceMismatch => 'X',
        };
        s.push(c);
    }
    s
}

/// Compute the aligned length on the reference using only Match and Deletion CIGAR operators.
pub fn get_aligned_length_md(cigar: &[(u32, CigarKind)]) -> i32 {
    cigar
        .iter()
        .filter(|&&(_, kind)| kind == CigarKind::Match || kind == CigarKind::Deletion)
        .map(|&(len, _)| len as i32)
        .sum()
}

/// Class for primary processing of SAM record and filtering reads that don't fit the requirements.
/// Ports RecordPreprocessor.java.
pub struct RecordPreprocessor {
    pub duplicates: HashSet<String>,
    pub first_matching_position: i32,
    pub total_reads: u32,
    pub duplicate_reads: u32,
}

impl RecordPreprocessor {
    pub fn new() -> Self {
        RecordPreprocessor {
            duplicates: HashSet::new(),
            first_matching_position: -1,
            total_reads: 0,
            duplicate_reads: 0,
        }
    }

    /// Preprocesses a record, returning true if it passes all filters.
    pub fn preprocess_record(
        &mut self,
        record: &SamRecord,
        config: &crate::config::Config,
    ) -> bool {
        // 1. Mapping quality filter
        if let Some(min_mapq) = config.map_qual {
            if let Some(mapq) = record.mapping_quality {
                if (mapq as u32) < min_mapq {
                    return false;
                }
            } else {
                // Missing MAPQ is treated as 255 which passes
                if 255 < min_mapq {
                    return false;
                }
            }
        }

        // 2. Secondary alignment filter (unless samfilter is explicitly 0)
        let is_secondary = (record.flags & 0x100) != 0;
        if is_secondary && config.flag_filter != 0 {
            return false;
        }

        // Unmapped reads filter
        if (record.flags & 0x4) != 0 {
            return false;
        }

        // 3. Sequence is not stored check
        if record.sequence.is_empty() || (record.sequence.len() == 1 && record.sequence[0] == b'*')
        {
            return false;
        }

        // 4. SAM flag filter
        if config.flag_filter != 0 && (record.flags & (config.flag_filter as u16)) != 0 {
            return false;
        }

        self.total_reads += 1;

        // 5. Duplicated reads filter (removeDuplicatedReads)
        if config.dedup {
            let alignment_start_val = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
            let mate_alignment_start_val =
                record.mate_alignment_start.map(|p| p.0 as i32).unwrap_or(0);

            if alignment_start_val != self.first_matching_position {
                self.duplicates.clear();
            }

            if mate_alignment_start_val < 10 {
                // POS-RNEXT-PNEXT
                let dup_key = format!(
                    "{}-{}-{}",
                    alignment_start_val,
                    get_mate_reference_name(record),
                    mate_alignment_start_val
                );
                if self.duplicates.contains(&dup_key) {
                    self.duplicate_reads += 1;
                    return false;
                }
                self.duplicates.insert(dup_key);
                self.first_matching_position = alignment_start_val;
            } else {
                let is_paired = (record.flags & 0x1) != 0;
                let is_mate_unmapped = (record.flags & 0x8) != 0;
                if is_paired && is_mate_unmapped {
                    // POS-CIGAR
                    let dup_key = format!(
                        "{}-{}",
                        alignment_start_val,
                        get_cigar_string(&record.cigar)
                    );
                    if self.duplicates.contains(&dup_key) {
                        self.duplicate_reads += 1;
                        return false;
                    }
                    self.duplicates.insert(dup_key);
                    self.first_matching_position = alignment_start_val;
                }
            }
        }

        true
    }
}

impl Default for RecordPreprocessor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_bam_open_and_query() {
        let bam_path = Path::new("/mnt/hbmon/VarDictRust/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam");
        let config = Config::default();

        let cached = open_reader(bam_path, &config).unwrap();
        let ref_seqs = cached.header.reference_sequences();
        assert!(
            !ref_seqs.is_empty(),
            "Reference sequences should not be empty"
        );

        let records =
            query_region(bam_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        assert!(!records.is_empty(), "Should load records from BAM file");
    }

    #[test]
    fn flag_filter() {
        let bam_path = Path::new("/mnt/hbmon/VarDictRust/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam");

        // 1. Without flag filter
        let config_no_filter = Config {
            flag_filter: 0,
            ..Default::default()
        };

        let records = query_region(
            bam_path,
            "chr7",
            Pos0(55259400),
            Pos0(55259600),
            &config_no_filter,
        )
        .unwrap();

        assert_eq!(
            records.len(),
            1327,
            "Total unfiltered record count must match samtools"
        );

        // 2. With flag filter 0x504
        let config_with_filter = Config {
            flag_filter: 0x504,
            ..Default::default()
        };

        let mut preprocessor = RecordPreprocessor::new();
        let mut passed_count = 0;

        for record in &records {
            if preprocessor.preprocess_record(record, &config_with_filter) {
                passed_count += 1;
            }
        }

        assert_eq!(
            passed_count, 1031,
            "Filtered count must match samtools view -c -F 0x504"
        );
    }

    #[test]
    fn test_cram_twin_parity() {
        // Detect samtools availability
        let mut has_samtools = std::process::Command::new("samtools")
            .arg("--version")
            .output()
            .is_ok();
        let mut use_conda = false;
        if !has_samtools {
            let conda_check = std::process::Command::new("conda")
                .args(["run", "-n", "vardict-oracle", "samtools", "--version"])
                .output();
            if conda_check.is_ok() {
                has_samtools = true;
                use_conda = true;
            }
        }
        if !has_samtools {
            println!("Skipping CRAM parity test: samtools not found");
            return;
        }

        let temp_dir = std::env::temp_dir().join("cram_test");
        let _ = std::fs::create_dir_all(&temp_dir);

        let ref_fa = temp_dir.join("ref.fa");
        let ref_fai = temp_dir.join("ref.fa.fai");
        let slice_sam = temp_dir.join("slice.sam");
        let new_sam_path = temp_dir.join("new_slice.sam");
        let slice_bam = temp_dir.join("slice.bam");
        let cram_path = temp_dir.join("test.cram");
        let bam_path = Path::new("/mnt/hbmon/VarDictRust/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam");

        // 1. Slice BAM to SAM
        let mut slice_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-h",
                "-o",
                slice_sam.to_str().unwrap(),
                bam_path.to_str().unwrap(),
                "chr7:55259401-55259600",
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-h",
                "-o",
                slice_sam.to_str().unwrap(),
                bam_path.to_str().unwrap(),
                "chr7:55259401-55259600",
            ]);
            c
        };
        assert!(slice_cmd.output().unwrap().status.success());

        // 2. Read slice.sam and write a new SAM file with a header containing only chr7
        let sam_content = std::fs::read_to_string(&slice_sam).unwrap();
        let mut new_sam = File::create(&new_sam_path).unwrap();
        let mut wrote_sq = false;
        for line in sam_content.lines() {
            if line.starts_with("@SQ") {
                if line.contains("SN:chr7") && !wrote_sq {
                    writeln!(new_sam, "{}", line).unwrap();
                    wrote_sq = true;
                }
            } else {
                writeln!(new_sam, "{}", line).unwrap();
            }
        }
        drop(new_sam);

        // 3. Convert SAM back to BAM (this re-encodes the reference IDs correctly)
        let mut convert_bam_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-b",
                "-o",
                slice_bam.to_str().unwrap(),
                new_sam_path.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-b",
                "-o",
                slice_bam.to_str().unwrap(),
                new_sam_path.to_str().unwrap(),
            ]);
            c
        };
        assert!(convert_bam_cmd.output().unwrap().status.success());

        // 4. Index the sliced BAM
        let mut slice_index = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "index",
                slice_bam.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["index", slice_bam.to_str().unwrap()]);
            c
        };
        assert!(slice_index.output().unwrap().status.success());

        // 5. Write a temporary mock reference FASTA with 159,138,663 N bases for chr7 to match BAM header chromosome length
        let mut fa_file = File::create(&ref_fa).unwrap();
        writeln!(fa_file, ">chr7").unwrap();
        let chunk = vec![b'N'; 100_000];
        for _ in 0..1591 {
            fa_file.write_all(&chunk).unwrap();
            fa_file.write_all(b"\n").unwrap();
        }
        let remainder = vec![b'N'; 38663];
        fa_file.write_all(&remainder).unwrap();
        fa_file.write_all(b"\n").unwrap();
        drop(fa_file);

        // 6. Run samtools faidx
        let mut faidx_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "faidx",
                ref_fa.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["faidx", ref_fa.to_str().unwrap()]);
            c
        };
        assert!(faidx_cmd.output().unwrap().status.success());

        // 7. Run samtools to convert slice_bam to CRAM
        let mut convert_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-C",
                "-T",
                ref_fa.to_str().unwrap(),
                "-o",
                cram_path.to_str().unwrap(),
                slice_bam.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-C",
                "-T",
                ref_fa.to_str().unwrap(),
                "-o",
                cram_path.to_str().unwrap(),
                slice_bam.to_str().unwrap(),
            ]);
            c
        };
        assert!(convert_cmd.output().unwrap().status.success());

        // 8. Run samtools index on CRAM file
        let mut index_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "index",
                cram_path.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["index", cram_path.to_str().unwrap()]);
            c
        };
        assert!(index_cmd.output().unwrap().status.success());

        // Now query both BAM and CRAM and verify parity
        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };

        // Query BAM
        let bam_records =
            query_region(&slice_bam, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        // Query CRAM
        let cram_records =
            query_region(&cram_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        assert_eq!(
            bam_records.len(),
            cram_records.len(),
            "BAM and CRAM must have the exact same number of reads"
        );

        for (i, (bam_rec, cram_rec)) in bam_records.iter().zip(cram_records.iter()).enumerate() {
            assert_eq!(bam_rec, cram_rec, "Parity mismatch at record {}", i);
        }

        // Cleanup
        let _ = std::fs::remove_file(&ref_fa);
        let _ = std::fs::remove_file(&ref_fai);
        let _ = std::fs::remove_file(&slice_sam);
        let _ = std::fs::remove_file(&new_sam_path);
        let _ = std::fs::remove_file(&slice_bam);
        let _ = std::fs::remove_file(slice_bam.with_extension("bam.bai"));
        let _ = std::fs::remove_file(&cram_path);
        let _ = std::fs::remove_file(cram_path.with_extension("cram.crai"));
    }
}
