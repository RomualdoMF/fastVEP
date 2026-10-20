//! Gene/transcript-keyed records for `.osa2`.
//!
//! A keyed record holds every (gene, transcript) entry of one (pos, ref, alt):
//! the shared JSON object, then each entry as `gene<US>transcript<US>{columns}`,
//! all joined by `<RS>` (control characters that JSON text never contains, and
//! not a newline, which the chunk's blob column uses to separate records).
//! Lookup compares only the two ID fields of each entry and parses JSON for the
//! matching entry alone.

use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Which keys a database carries; stored under `keys` in `metadata.json`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KeyInfo {
    pub gene: bool,
    pub transcript: bool,
    /// `ensembl` or `refseq`: the ID namespace of the stored keys.
    pub key_type: String,
    /// Annotation columns of the source, in input order. `annotate` lists them
    /// in the VCF CSQ header; empty for databases built before this existed.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub columns: Vec<String>,
}

/// One entry of a keyed record, as the builder collects it.
pub struct KeyedEntry {
    pub gene: Option<String>,
    pub transcript: Option<String>,
    /// JSON object text holding the per-entry columns.
    pub columns: String,
}

/// Drop an Ensembl/RefSeq version suffix (`ENST000001.5` -> `ENST000001`).
pub fn strip_version(id: &str) -> &str {
    match id.rsplit_once('.') {
        Some((base, v)) if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) => base,
        _ => id,
    }
}

fn same_id(a: &str, b: &str) -> bool {
    strip_version(a).eq_ignore_ascii_case(strip_version(b))
}

/// Whether two keys name the same thing under `key_type`.
fn same_key(key_type: &str, a: &str, b: &str) -> bool {
    if key_type == "symbol" {
        a.eq_ignore_ascii_case(b)
    } else {
        same_id(a, b)
    }
}

/// True when `id` looks like it belongs to `key_type`'s namespace. `symbol`
/// accepts any non-empty name and is compared exactly (no version suffix:
/// `AL627309.1` and `AL627309.2` are different genes).
pub fn id_matches_type(id: &str, is_transcript: bool, key_type: &str) -> bool {
    let ensembl = id.starts_with("ENS");
    match key_type {
        // Gene symbols: any name will do.
        "symbol" => !id.is_empty(),
        "ensembl" => ensembl,
        "refseq" if is_transcript => ["NM_", "NR_", "XM_", "XR_"]
            .iter()
            .any(|p| id.starts_with(p)),
        "refseq" => !ensembl,
        _ => false,
    }
}

const ENTRY_SEP: char = '\u{1e}';
const FIELD_SEP: char = '\u{1f}';

/// IDs sit between separators, so control characters are blanked.
fn clean_id(id: &str) -> String {
    id.chars().map(|c| if c.is_control() { ' ' } else { c }).collect()
}

/// Serialize a keyed record; `shared` is a JSON object text (single line).
pub fn build_blob(shared: &str, entries: &[KeyedEntry]) -> String {
    let mut out = String::with_capacity(shared.len() + entries.len() * 256);
    out.push_str(shared);
    for e in entries {
        out.push(ENTRY_SEP);
        out.push_str(&clean_id(e.gene.as_deref().unwrap_or("")));
        out.push(FIELD_SEP);
        out.push_str(&clean_id(e.transcript.as_deref().unwrap_or("")));
        out.push(FIELD_SEP);
        out.push_str(&e.columns);
    }
    out
}

/// Pick the entry matching the query and return `shared` merged with it.
///
/// With `query == None` (an unkeyed call) the first entry is returned. With a
/// query, every key the database carries must match; no match gives `None`.
/// Entries that do not match are never parsed.
pub fn select(
    blob: &str,
    keys: &KeyInfo,
    query: Option<(Option<&str>, Option<&str>)>,
) -> Option<String> {
    let mut lines = blob.split(ENTRY_SEP);
    let shared = lines.next()?;
    let columns = lines.find_map(|line| {
        let (gene, rest) = line.split_once(FIELD_SEP)?;
        let (tx, columns) = rest.split_once(FIELD_SEP)?;
        let matches = match query {
            None => true,
            Some((want_gene, want_tx)) => {
                let ok = |enabled: bool, have: &str, want: Option<&str>| {
                    !enabled || want.is_some_and(|w| !have.is_empty() && same_key(&keys.key_type, have, w))
                };
                ok(keys.gene, gene, want_gene) && ok(keys.transcript, tx, want_tx)
            }
        };
        matches.then_some(columns)
    })?;
    // Nothing to merge when either side is empty.
    if shared == "{}" || shared.is_empty() {
        return Some(columns.to_string());
    }
    if columns == "{}" {
        return Some(shared.to_string());
    }
    let (Value::Object(mut merged), Value::Object(cols)) = (
        serde_json::from_str::<Value>(shared).ok()?,
        serde_json::from_str::<Value>(columns).ok()?,
    ) else {
        return None;
    };
    merged.extend(cols);
    serde_json::to_string(&Value::Object(merged)).ok()
}

/// One input line's worth of entries for a (pos, ref, alt).
pub struct KeyedRow {
    pub chrom_idx: u16,
    pub position: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    /// JSON object text of the columns common to all entries.
    pub shared: String,
    pub entries: Vec<KeyedEntry>,
}

struct Group {
    ref_allele: String,
    alt_allele: String,
    shared: String,
    entries: Vec<KeyedEntry>,
    seen: std::collections::HashSet<(String, String)>,
}

/// Folds position-sorted [`KeyedRow`]s into one record per (pos, ref, alt).
///
/// Entries with a repeated (gene, transcript) keep the first; `shared` comes
/// from the first row of the group.
pub struct Grouper<I> {
    rows: I,
    /// Drop `.N` suffixes before telling two IDs apart (not for gene symbols).
    strip_versions: bool,
    at: Option<(u16, u32)>,
    groups: Vec<Group>,
    ready: std::collections::VecDeque<crate::common::AnnotationRecord>,
    done: bool,
}

impl<I> Grouper<I> {
    pub fn new(rows: I) -> Self {
        Self::with_versions(rows, true)
    }

    /// As [`new`](Self::new) for gene symbols, which are compared exactly.
    pub fn exact(rows: I) -> Self {
        Self::with_versions(rows, false)
    }

    fn with_versions(rows: I, strip_versions: bool) -> Self {
        Self {
            rows,
            strip_versions,
            at: None,
            groups: Vec::new(),
            ready: Default::default(),
            done: false,
        }
    }

    fn flush(&mut self) {
        let Some((chrom_idx, position)) = self.at else { return };
        for g in self.groups.drain(..) {
            self.ready.push_back(crate::common::AnnotationRecord {
                chrom_idx,
                position,
                ref_allele: g.ref_allele,
                alt_allele: g.alt_allele,
                json: build_blob(&g.shared, &g.entries),
            });
        }
    }

    fn add(&mut self, row: KeyedRow) {
        if self.at != Some((row.chrom_idx, row.position)) {
            self.flush();
            self.at = Some((row.chrom_idx, row.position));
        }
        let gi = match self
            .groups
            .iter()
            .position(|g| g.ref_allele == row.ref_allele && g.alt_allele == row.alt_allele)
        {
            Some(i) => i,
            None => {
                self.groups.push(Group {
                    ref_allele: row.ref_allele,
                    alt_allele: row.alt_allele,
                    shared: row.shared,
                    entries: Vec::new(),
                    seen: Default::default(),
                });
                self.groups.len() - 1
            }
        };
        let group = &mut self.groups[gi];
        for e in row.entries {
            let strip = self.strip_versions;
            let norm = |v: &Option<String>| {
                v.as_deref().map_or(String::new(), |v| {
                    (if strip { strip_version(v) } else { v }).to_ascii_lowercase()
                })
            };
            if group.seen.insert((norm(&e.gene), norm(&e.transcript))) {
                group.entries.push(e);
            }
        }
    }
}

impl<I: Iterator<Item = anyhow::Result<KeyedRow>>> Iterator for Grouper<I> {
    type Item = anyhow::Result<crate::common::AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(rec) = self.ready.pop_front() {
                return Some(Ok(rec));
            }
            if self.done {
                return None;
            }
            match self.rows.next() {
                Some(Ok(row)) => self.add(row),
                Some(Err(e)) => {
                    self.done = true;
                    return Some(Err(e));
                }
                None => {
                    self.done = true;
                    self.flush();
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn keys(gene: bool, transcript: bool) -> KeyInfo {
        KeyInfo {
            gene,
            transcript,
            key_type: "ensembl".into(),
            columns: Vec::new(),
        }
    }

    fn blob() -> String {
        build_blob(
            r#"{"CADD_phred":"12"}"#,
            &[
                KeyedEntry {
                    gene: Some("ENSG1".into()),
                    transcript: Some("ENST1".into()),
                    columns: r#"{"SIFT_score":"0.1"}"#.into(),
                },
                KeyedEntry {
                    gene: Some("ENSG1".into()),
                    transcript: Some("ENST2".into()),
                    columns: r#"{"SIFT_score":"0.9"}"#.into(),
                },
            ],
        )
    }

    #[test]
    fn versions_are_ignored_on_both_sides() {
        assert_eq!(strip_version("ENST0001.12"), "ENST0001");
        assert_eq!(strip_version("NM_000546.6"), "NM_000546");
        assert_eq!(strip_version("ENST0001"), "ENST0001");
        assert!(same_id("ENST2.3", "enst2"));
    }

    #[test]
    fn transcript_key_picks_its_entry_and_merges_shared() {
        let j = select(&blob(), &keys(false, true), Some((None, Some("ENST2.7")))).unwrap();
        assert!(j.contains("\"SIFT_score\":\"0.9\"") && j.contains("CADD_phred"), "{j}");
    }

    #[test]
    fn no_match_is_empty() {
        assert!(select(&blob(), &keys(false, true), Some((None, Some("ENST9")))).is_none());
        // A key the database needs but the query lacks cannot match.
        assert!(select(&blob(), &keys(true, true), Some((None, Some("ENST1")))).is_none());
    }

    #[test]
    fn gene_only_key_keeps_the_first_entry_of_the_gene() {
        let j = select(&blob(), &keys(true, false), Some((Some("ENSG1"), None))).unwrap();
        assert!(j.contains("\"SIFT_score\":\"0.1\""), "{j}");
    }

    #[test]
    fn non_matching_entries_are_never_parsed() {
        // The first entry's payload is not even valid JSON; a lookup for the
        // second must still succeed, which it can only do by skipping it.
        let blob = "{\"s\":\"1\"}\u{1e}ENSG1\u{1f}ENST1\u{1f}{not json\u{1e}ENSG1\u{1f}ENST2\u{1f}{\"v\":\"2\"}";
        let j = select(blob, &keys(true, true), Some((Some("ENSG1"), Some("ENST2")))).unwrap();
        assert!(j.contains("\"v\":\"2\"") && j.contains("\"s\":\"1\""), "{j}");
        assert!(select(blob, &keys(true, true), Some((Some("ENSG1"), Some("ENST1")))).is_none());
    }

    #[test]
    fn unkeyed_call_returns_the_first_entry() {
        let j = select(&blob(), &keys(true, true), None).unwrap();
        assert!(j.contains("\"SIFT_score\":\"0.1\""), "{j}");
    }

    #[test]
    fn grouper_merges_rows_and_keeps_the_first_duplicate() {
        let row = |t: &str, v: &str, shared: &str| {
            Ok(KeyedRow {
                chrom_idx: 0,
                position: 10,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                shared: shared.into(),
                entries: vec![KeyedEntry {
                    gene: None,
                    transcript: Some(t.into()),
                    columns: format!("{{\"v\":\"{v}\"}}"),
                }],
            })
        };
        let rows = vec![row("T1", "a", "{\"s\":\"1\"}"), row("T2.3", "b", "{}"), row("t2", "c", "{}")];
        let recs: Vec<_> = Grouper::new(rows.into_iter()).collect::<anyhow::Result<_>>().unwrap();
        assert_eq!(recs.len(), 1);
        let k = keys(false, true);
        let hit = |t| select(&recs[0].json, &k, Some((None, Some(t)))).unwrap();
        assert!(hit("T2").contains("\"v\":\"b\""), "first T2 wins");
        assert!(hit("T1").contains("\"s\":\"1\""), "shared from the first row");
    }

    #[test]
    fn gene_symbols_are_compared_exactly() {
        let symbols = KeyInfo { gene: true, transcript: false, key_type: "symbol".into(), columns: vec![] };
        let blob = build_blob(
            "{}",
            &[
                KeyedEntry { gene: Some("AL627309.1".into()), transcript: None, columns: "{\"v\":\"1\"}".into() },
                KeyedEntry { gene: Some("AL627309.2".into()), transcript: None, columns: "{\"v\":\"2\"}".into() },
            ],
        );
        // Under a version-stripping key type both would be the same gene.
        let two = select(&blob, &symbols, Some((Some("AL627309.2"), None))).unwrap();
        assert!(two.contains("\"v\":\"2\""), "{two}");
        assert!(select(&blob, &symbols, Some((Some("al627309.1"), None))).unwrap().contains("\"v\":\"1\""));
        assert!(select(&blob, &symbols, Some((Some("AL627309"), None))).is_none());
        // The grouper keeps both, where ensembl-style IDs would have merged them.
        let row = |g: &str| Ok(KeyedRow {
            chrom_idx: 0, position: 1, ref_allele: "A".into(), alt_allele: "G".into(), shared: "{}".into(),
            entries: vec![KeyedEntry { gene: Some(g.into()), transcript: None, columns: "{}".into() }],
        });
        let recs: Vec<_> = Grouper::exact(vec![row("AL627309.1"), row("AL627309.2")].into_iter())
            .collect::<anyhow::Result<_>>().unwrap();
        assert_eq!(recs[0].json.matches('\u{1e}').count(), 2);
        let merged: Vec<_> = Grouper::new(vec![row("AL627309.1"), row("AL627309.2")].into_iter())
            .collect::<anyhow::Result<_>>().unwrap();
        assert_eq!(merged[0].json.matches('\u{1e}').count(), 1);
        assert!(id_matches_type("AL627309.1", false, "symbol") && !id_matches_type("", false, "symbol"));
    }

    #[test]
    fn id_namespaces() {
        assert!(id_matches_type("ENST1", true, "ensembl"));
        assert!(!id_matches_type("ENST1", true, "refseq"));
        assert!(id_matches_type("NM_000546.6", true, "refseq"));
        assert!(!id_matches_type("NM_1", true, "ensembl"));
    }
}
