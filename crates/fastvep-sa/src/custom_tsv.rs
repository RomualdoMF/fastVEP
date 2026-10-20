//! `custom_tsv`: user-supplied tab-separated annotations for `.osa2`.
//!
//! The file has a header row naming a chromosome, position, ref and alt column;
//! every other column becomes a JSON key. Optional gene / transcript columns
//! turn the database into a keyed one (see [`crate::keyed`]).

use crate::common::AnnotationRecord;
use crate::keyed::{id_matches_type, KeyInfo, KeyedEntry, KeyedRow};
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, Lines};

const CHR: &[&str] = &["chr", "chrom", "chromosome"];
const POS: &[&str] = &["pos", "position", "pos(1-based)"];
const REF: &[&str] = &["ref", "reference", "ref_allele"];
const ALT: &[&str] = &["alt", "alternate", "alt_allele"];
/// Separator of the values inside one cell (`ENST1;ENST2`).
const LIST_SEP: char = ';';

/// Metadata for a custom TSV database named `name`.
pub fn custom_tsv_osa2_metadata(name: &str, assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: name.into(),
        version: "user-supplied".into(),
        assembly: assembly.into(),
        json_key: name.into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("Custom TSV annotations for {assembly}"),
    }
}

/// Columns holding the gene / transcript ID, and the ID namespace.
pub struct TsvKeys {
    pub transcript: Option<String>,
    pub gene: Option<String>,
    /// `ensembl` or `refseq`.
    pub key_type: String,
}

impl TsvKeys {
    /// The `keys` block stored in `metadata.json`; `columns` are the stored
    /// annotation columns (see [`custom_tsv_columns`]).
    pub fn info(&self, columns: Vec<String>) -> KeyInfo {
        KeyInfo {
            gene: self.gene.is_some(),
            transcript: self.transcript.is_some(),
            key_type: self.key_type.clone(),
            columns,
        }
    }
}

struct Layout {
    chr: usize,
    pos: usize,
    ref_col: usize,
    alt: usize,
    names: Vec<String>,
    /// Columns stored as annotation values.
    info: Vec<usize>,
    tx: Option<usize>,
    gene: Option<usize>,
}

impl Layout {
    fn from_header(header: &str, info_fields: &[String], keys: Option<&TsvKeys>) -> Result<Self> {
        let names: Vec<String> = header
            .split('\t')
            .map(|h| h.trim().trim_start_matches('#').to_string())
            .collect();
        let find = |aliases: &[&str], what: &str| -> Result<usize> {
            names
                .iter()
                .position(|n| aliases.iter().any(|a| n.eq_ignore_ascii_case(a)))
                .with_context(|| format!("custom TSV header has no {what} column ({aliases:?})"))
        };
        let (chr, pos) = (find(CHR, "chromosome")?, find(POS, "position")?);
        let (ref_col, alt) = (find(REF, "ref")?, find(ALT, "alt")?);
        let by_name = |name: &str| -> Result<usize> {
            names
                .iter()
                .position(|n| n == name)
                .with_context(|| format!("column '{name}' not found in the TSV header"))
        };
        let info = if info_fields.is_empty() {
            (0..names.len())
                .filter(|i| ![chr, pos, ref_col, alt].contains(i))
                .collect()
        } else {
            info_fields.iter().map(|f| by_name(f)).collect::<Result<_>>()?
        };
        let key_col = |c: Option<&String>| c.map(|c| by_name(c)).transpose();
        Ok(Self {
            tx: key_col(keys.and_then(|k| k.transcript.as_ref()))?,
            gene: key_col(keys.and_then(|k| k.gene.as_ref()))?,
            chr,
            pos,
            ref_col,
            alt,
            names,
            info,
        })
    }
}

/// A data row whose chromosome and position already parsed.
struct Row {
    chrom_idx: u16,
    position: u32,
    line: String,
}

/// Reads the header, then yields data rows of known chromosomes.
struct Source<'a, R: BufRead> {
    lines: Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    info_fields: &'a [String],
    keys: Option<&'a TsvKeys>,
    layout: Option<Layout>,
}

impl<'a, R: BufRead> Source<'a, R> {
    fn new(
        reader: R,
        chrom_to_idx: &'a HashMap<String, u16>,
        info_fields: &'a [String],
        keys: Option<&'a TsvKeys>,
    ) -> Self {
        Self {
            lines: reader.lines(),
            chrom_to_idx,
            info_fields,
            keys,
            layout: None,
        }
    }

    fn next_row(&mut self) -> Option<Result<Row>> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(anyhow::Error::new(e).context("Reading custom TSV"))),
            };
            let line = line.trim_end_matches('\r');
            if line.is_empty() || line.starts_with("##") {
                continue;
            }
            // The header is the first non-`##` line.
            let Some(layout) = &self.layout else {
                match Layout::from_header(line, self.info_fields, self.keys) {
                    Ok(l) => self.layout = Some(l),
                    Err(e) => return Some(Err(e)),
                }
                continue;
            };
            let mut fields = line.split('\t');
            let chrom = fields.nth(layout.chr)?;
            let chrom = if chrom.starts_with("chr") {
                chrom.to_string()
            } else {
                format!("chr{chrom}")
            };
            let Some(&chrom_idx) = self.chrom_to_idx.get(&chrom) else {
                continue;
            };
            let Some(Ok(position)) = line.split('\t').nth(layout.pos).map(str::parse::<u32>)
            else {
                continue;
            };
            return Some(Ok(Row {
                chrom_idx,
                position,
                line: line.to_string(),
            }));
        }
    }
}

/// Names of the annotation columns a build stores, read from the header alone.
pub fn custom_tsv_columns<R: BufRead>(reader: R, info_fields: &[String]) -> Result<Vec<String>> {
    for line in reader.lines() {
        let line = line.context("Reading custom TSV")?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() || line.starts_with("##") {
            continue;
        }
        let layout = Layout::from_header(line, info_fields, None)?;
        return Ok(layout.info.iter().map(|&c| layout.names[c].clone()).collect());
    }
    anyhow::bail!("custom TSV has no header line")
}

fn is_missing(v: &str) -> bool {
    v.is_empty() || v == "."
}

/// Append `"key":"value"`, escaping only when a character needs it.
fn push_kv(out: &mut String, key: &str, value: &str) {
    for (i, part) in [key, value].into_iter().enumerate() {
        out.push('"');
        if part.bytes().all(|b| b >= 0x20 && b != b'"' && b != b'\\') {
            out.push_str(part);
        } else {
            out.push_str(&crate::common::escape_json(part));
        }
        out.push_str(if i == 0 { "\":" } else { "\"" });
    }
}

fn push_sep(out: &mut String) {
    out.push(if out.is_empty() { '{' } else { ',' });
}

fn finish_object(mut out: String) -> String {
    out.push_str(if out.is_empty() { "{}" } else { "}" });
    out
}

/// True when the rows can stream into the v2 writer without a sort: each
/// chromosome appears in one block, with positions not decreasing inside it.
/// The chromosome order itself does not matter (dbNSFP lists 1, 10, 11, ...).
pub fn custom_tsv_is_sorted<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<bool> {
    let mut src = Source::new(reader, chrom_to_idx, &[], None);
    let mut done = std::collections::HashSet::new();
    let (mut chrom, mut last) = (None, 0u32);
    while let Some(row) = src.next_row() {
        let row = row?;
        if chrom != Some(row.chrom_idx) {
            if !done.insert(row.chrom_idx) {
                return Ok(false);
            }
            (chrom, last) = (Some(row.chrom_idx), 0);
        }
        if row.position < last {
            return Ok(false);
        }
        last = row.position;
    }
    Ok(true)
}

/// Stream a TSV as one record per row, holding every selected column.
///
/// Missing values (`.` or empty) are left out. Holds one row at a time.
pub fn iter_custom_tsv<'a, R: BufRead + 'a>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    name: &'a str,
    info_fields: &'a [String],
) -> impl Iterator<Item = Result<AnnotationRecord>> + 'a {
    let mut src = Source::new(reader, chrom_to_idx, info_fields, None);
    std::iter::from_fn(move || loop {
        let row = match src.next_row()? {
            Ok(r) => r,
            Err(e) => return Some(Err(e)),
        };
        let layout = src.layout.as_ref()?;
        let fields: Vec<&str> = row.line.split('\t').collect();
        let (Some(&ref_allele), Some(&alt_allele)) =
            (fields.get(layout.ref_col), fields.get(layout.alt))
        else {
            continue;
        };
        if is_missing(alt_allele) || alt_allele == "*" {
            continue;
        }
        let mut json = String::with_capacity(row.line.len() + 32);
        for &c in &layout.info {
            let Some(v) = fields.get(c).filter(|v| !is_missing(v)) else {
                continue;
            };
            push_sep(&mut json);
            push_kv(&mut json, &layout.names[c], v);
        }
        if json.is_empty() {
            push_sep(&mut json);
            push_kv(&mut json, "source", name);
        }
        return Some(Ok(AnnotationRecord {
            chrom_idx: row.chrom_idx,
            position: row.position,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            json: finish_object(json),
        }));
    })
}

/// Keyed form of [`iter_custom_tsv`]: one [`KeyedRow`] per input row.
///
/// A key cell with N IDs (joined by `;`) makes N entries. Any other selected
/// cell with exactly N values is split the same way; the rest are shared.
/// Rows without the IDs, or with differing gene / transcript counts, are skipped.
pub fn iter_custom_tsv_keyed<'a, R: BufRead + 'a>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    info_fields: &'a [String],
    keys: &'a TsvKeys,
) -> impl Iterator<Item = Result<KeyedRow>> + 'a {
    let mut src = Source::new(reader, chrom_to_idx, info_fields, Some(keys));
    let mut warned = false;
    std::iter::from_fn(move || loop {
        let row = match src.next_row()? {
            Ok(r) => r,
            Err(e) => return Some(Err(e)),
        };
        let layout = src.layout.as_ref()?;
        let fields: Vec<&str> = row.line.split('\t').collect();
        let (Some(&ref_allele), Some(&alt_allele)) =
            (fields.get(layout.ref_col), fields.get(layout.alt))
        else {
            continue;
        };
        let ids = |c: Option<usize>| -> Option<Vec<&str>> {
            c.and_then(|c| fields.get(c))
                .map(|v| v.split(LIST_SEP).map(str::trim).collect())
        };
        let (tx, gene) = (ids(layout.tx), ids(layout.gene));
        let n = match (&tx, &gene) {
            (Some(t), Some(g)) if t.len() != g.len() => 0,
            (Some(l), _) | (_, Some(l)) => l.len(),
            _ => 0,
        };
        if n == 0 || is_missing(alt_allele) || alt_allele == "*" {
            if !warned && n == 0 {
                eprintln!("warning: skipping custom TSV rows without matching key values");
                warned = true;
            }
            continue;
        }

        // Selected cells split once; a cell is per-entry when it has n values.
        let split: Vec<(&str, Vec<&str>)> = layout
            .info
            .iter()
            .filter_map(|&c| fields.get(c).map(|v| (layout.names[c].as_str(), v.split(LIST_SEP).collect())))
            .collect();
        let mut shared = String::new();
        for (k, parts) in split.iter().filter(|(_, p)| p.len() != n) {
            let whole = parts.join(&LIST_SEP.to_string());
            if !is_missing(&whole) {
                push_sep(&mut shared);
                push_kv(&mut shared, k, &whole);
            }
        }

        let mut entries = Vec::with_capacity(n);
        for i in 0..n {
            let (t, g) = (tx.as_ref().map(|l| l[i]), gene.as_ref().map(|l| l[i]));
            if [t, g].into_iter().flatten().any(is_missing) {
                continue;
            }
            for (id, is_tx) in [(t, true), (g, false)] {
                if let Some(id) = id.filter(|id| !id_matches_type(id, is_tx, &keys.key_type)) {
                    return Some(Err(anyhow::anyhow!(
                        "{} '{}' is not a {} ID (--key-type {})",
                        if is_tx { "transcript" } else { "gene" },
                        id,
                        keys.key_type,
                        keys.key_type
                    )));
                }
            }
            let mut columns = String::new();
            for (k, parts) in split.iter().filter(|(_, p)| p.len() == n) {
                if !is_missing(parts[i]) {
                    push_sep(&mut columns);
                    push_kv(&mut columns, k, parts[i]);
                }
            }
            entries.push(KeyedEntry {
                gene: g.map(str::to_string),
                transcript: t.map(str::to_string),
                columns: finish_object(columns),
            });
        }
        return Some(Ok(KeyedRow {
            chrom_idx: row.chrom_idx,
            position: row.position,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            shared: finish_object(shared),
            entries,
        }));
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chroms() -> HashMap<String, u16> {
        HashMap::from([("chr1".to_string(), 0u16), ("chr2".to_string(), 1)])
    }

    const TSV: &str = "#chr\tpos(1-based)\tref\talt\tscore\tgene\ttx\tsift\n\
        1\t100\tA\tG\t.\tENSG1;ENSG1\tENST1;ENST2\t0.1;0.9\n\
        2\t50\tC\tT\t7\tENSG2\tENST3\t0.5\n";

    #[test]
    fn plain_rows_keep_every_non_missing_column() {
        let recs: Vec<_> = iter_custom_tsv(TSV.as_bytes(), &chroms(), "t", &[])
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].chrom_idx, 0);
        assert!(!recs[0].json.contains("score"), "{}", recs[0].json);
        assert!(recs[1].json.contains("\"score\":\"7\""), "{}", recs[1].json);
    }

    #[test]
    fn info_fields_select_columns() {
        let f = vec!["sift".to_string()];
        let recs: Vec<_> = iter_custom_tsv(TSV.as_bytes(), &chroms(), "t", &f)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs[1].json, "{\"sift\":\"0.5\"}");
    }

    #[test]
    fn columns_come_from_the_header() {
        let all = custom_tsv_columns(TSV.as_bytes(), &[]).unwrap();
        assert_eq!(all, ["score", "gene", "tx", "sift"]);
        let some = custom_tsv_columns(TSV.as_bytes(), &["sift".to_string()]).unwrap();
        assert_eq!(some, ["sift"]);
    }

    #[test]
    fn missing_position_column_is_named() {
        let bad = "chr\tref\talt\n1\tA\tG\n";
        let err = iter_custom_tsv(bad.as_bytes(), &chroms(), "t", &[])
            .collect::<Result<Vec<_>>>()
            .unwrap_err();
        assert!(err.to_string().contains("position"), "{err}");
    }

    #[test]
    fn keyed_rows_split_parallel_cells() {
        let keys = TsvKeys {
            transcript: Some("tx".into()),
            gene: Some("gene".into()),
            key_type: "ensembl".into(),
        };
        assert_eq!(keys.info(vec!["a".into()]).columns, ["a"]);
        let rows: Vec<_> = iter_custom_tsv_keyed(TSV.as_bytes(), &chroms(), &[], &keys)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(rows[0].entries.len(), 2);
        assert_eq!(rows[0].entries[1].transcript.as_deref(), Some("ENST2"));
        assert!(rows[0].entries[1].columns.contains("\"sift\":\"0.9\""));
        assert_eq!(rows[1].entries.len(), 1);
    }

    #[test]
    fn sortedness_is_detected() {
        assert!(custom_tsv_is_sorted(TSV.as_bytes(), &chroms()).unwrap());
        // Any chromosome order streams, as long as each is one block.
        let reversed = "chr\tpos\tref\talt\n2\t5\tA\tG\n2\t9\tA\tG\n1\t9\tA\tG\n";
        assert!(custom_tsv_is_sorted(reversed.as_bytes(), &chroms()).unwrap());
        // A chromosome that returns, or a position that goes back, does not.
        let reopened = "chr\tpos\tref\talt\n1\t5\tA\tG\n2\t9\tA\tG\n1\t9\tA\tG\n";
        assert!(!custom_tsv_is_sorted(reopened.as_bytes(), &chroms()).unwrap());
        let backwards = "chr\tpos\tref\talt\n1\t9\tA\tG\n1\t5\tA\tG\n";
        assert!(!custom_tsv_is_sorted(backwards.as_bytes(), &chroms()).unwrap());
    }
}
