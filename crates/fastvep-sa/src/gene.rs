//! Gene-level annotation reader/writer (.oga files).
//!
//! Used for gene-keyed databases (OMIM, gnomAD gene scores, ClinGen).
//! Annotations are looked up by gene symbol.

use crate::common::{GeneRecord, MAX_INDEX_PAYLOAD, OGA_MAGIC, OGA_TRANSCRIPT_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use fastvep_cache::annotation::GeneAnnotationProvider;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};

/// Header for .oga files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneHeader {
    pub schema_version: u16,
    pub json_key: String,
    pub name: String,
    pub version: String,
    pub assembly: String,
}

/// In-memory gene annotation database loaded from an .oga file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneIndex {
    pub header: GeneHeader,
    /// Gene symbol -> list of JSON annotations (transcript ID when `transcript_keyed`).
    pub genes: HashMap<String, Vec<String>>,
    /// Records are filed under transcript IDs. Not part of the payload: it is
    /// told by the file's magic, so symbol-keyed files keep their old layout.
    #[serde(skip)]
    pub transcript_keyed: bool,
    /// Record keys of a transcript-keyed file, in column order.
    #[serde(skip)]
    pub columns: Vec<String>,
}

/// Payload of a transcript-keyed file: the index plus its column list.
#[derive(Serialize)]
struct TranscriptPayloadRef<'a> {
    index: &'a GeneIndex,
    columns: &'a [String],
}

#[derive(Deserialize)]
struct TranscriptPayload {
    index: GeneIndex,
    columns: Vec<String>,
}

/// An ID without its `.version` suffix (`ENST00000318602.11`, `NM_130786.4`).
pub fn strip_version(id: &str) -> &str {
    id.rsplit_once('.').filter(|(_, v)| v.chars().all(|c| c.is_ascii_digit()) && !v.is_empty()).map_or(id, |(b, _)| b)
}

impl GeneIndex {
    /// Create a new empty gene index.
    pub fn new(header: GeneHeader) -> Self {
        Self {
            header,
            genes: HashMap::new(),
            transcript_keyed: false,
            columns: Vec::new(),
        }
    }

    /// A new index whose records are filed under transcript IDs.
    pub fn new_transcript_keyed(header: GeneHeader, columns: Vec<String>) -> Self {
        Self { transcript_keyed: true, columns, ..Self::new(header) }
    }

    /// Add a gene record.
    pub fn add(&mut self, record: GeneRecord) {
        self.genes
            .entry(record.gene_symbol)
            .or_default()
            .push(record.json);
    }

    /// Look up annotations for a gene symbol.
    pub fn get(&self, gene_symbol: &str) -> Option<&Vec<String>> {
        self.genes.get(gene_symbol)
    }

    /// Number of genes in the index.
    pub fn gene_count(&self) -> usize {
        self.genes.len()
    }

    /// Serialize to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(if self.transcript_keyed { OGA_TRANSCRIPT_MAGIC } else { OGA_MAGIC })?;
        writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        let data = if self.transcript_keyed {
            bincode::serialize(&TranscriptPayloadRef { index: self, columns: &self.columns })?
        } else {
            bincode::serialize(self)?
        };
        writer.write_all(&(data.len() as u64).to_le_bytes())?;
        writer.write_all(&data)?;
        Ok(())
    }

    /// Deserialize from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        let transcript_keyed = &magic == OGA_TRANSCRIPT_MAGIC;
        if &magic != OGA_MAGIC && !transcript_keyed {
            anyhow::bail!("Invalid OGA magic");
        }
        let mut ver = [0u8; 2];
        reader.read_exact(&mut ver)?;
        if u16::from_le_bytes(ver) != SCHEMA_VERSION {
            anyhow::bail!("Unsupported OGA schema version");
        }
        let mut len_bytes = [0u8; 8];
        reader.read_exact(&mut len_bytes)?;
        let len_u64 = u64::from_le_bytes(len_bytes);
        if len_u64 > MAX_INDEX_PAYLOAD {
            anyhow::bail!(
                "OGA payload size {} exceeds limit {}",
                len_u64,
                MAX_INDEX_PAYLOAD
            );
        }
        let len: usize = len_u64
            .try_into()
            .map_err(|_| anyhow::anyhow!("OGA payload size {} exceeds usize", len_u64))?;
        let mut data = vec![0u8; len];
        reader.read_exact(&mut data)?;
        if transcript_keyed {
            let TranscriptPayload { mut index, columns } = bincode::deserialize(&data)?;
            index.transcript_keyed = true;
            index.columns = columns;
            return Ok(index);
        }
        let index: GeneIndex = bincode::deserialize(&data)?;
        Ok(index)
    }
}

impl GeneAnnotationProvider for GeneIndex {
    fn name(&self) -> &str {
        &self.header.name
    }

    fn json_key(&self) -> &str {
        &self.header.json_key
    }

    fn annotate_gene(&self, gene_symbol: &str) -> Result<Option<String>> {
        match self.genes.get(gene_symbol) {
            Some(jsons) if jsons.len() == 1 => Ok(Some(jsons[0].clone())),
            Some(jsons) => {
                // Multiple annotations: return as JSON array
                let array = format!("[{}]", jsons.join(","));
                Ok(Some(array))
            }
            None => Ok(None),
        }
    }

    fn is_transcript_keyed(&self) -> bool {
        self.transcript_keyed
    }

    fn annotate_transcript(&self, transcript_id: &str) -> Result<Option<String>> {
        // The first record of a transcript wins; a duplicate row is not a second answer.
        Ok(self.genes.get(strip_version(transcript_id)).and_then(|v| v.first().cloned()))
    }

    fn csq_columns(&self) -> &[String] {
        &self.columns
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_round_trip() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "omim".into(),
            name: "OMIM".into(),
            version: "2024-01".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = GeneIndex::new(header);
        index.add(GeneRecord {
            gene_symbol: "BRCA1".into(),
            json: r#"{"mim":113705}"#.into(),
        });
        index.add(GeneRecord {
            gene_symbol: "TP53".into(),
            json: r#"{"mim":191170}"#.into(),
        });

        assert_eq!(index.gene_count(), 2);

        // Serialize and deserialize
        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        let loaded = GeneIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();

        assert_eq!(loaded.header.json_key, "omim");
        assert_eq!(loaded.gene_count(), 2);
        assert_eq!(loaded.get("BRCA1").unwrap().len(), 1);
    }

    #[test]
    fn test_gene_provider() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "omim".into(),
            name: "OMIM".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = GeneIndex::new(header);
        index.add(GeneRecord {
            gene_symbol: "BRCA1".into(),
            json: r#"{"mim":113705}"#.into(),
        });

        let result = index.annotate_gene("BRCA1").unwrap();
        assert!(result.is_some());
        assert!(result.unwrap().contains("113705"));

        let result = index.annotate_gene("NONEXISTENT").unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn transcript_keyed_file_round_trips_and_ignores_versions() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "constraint".into(),
            name: "constraint".into(),
            version: "1".into(),
            assembly: "GRCh38".into(),
        };
        let mut index = GeneIndex::new_transcript_keyed(header, vec!["pLI".into()]);
        index.add(GeneRecord { gene_symbol: "ENST00000318602".into(), json: r#"{"pLI":"1"}"#.into() });
        index.add(GeneRecord { gene_symbol: "NM_130786".into(), json: r#"{"pLI":"2"}"#.into() });

        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        let loaded = GeneIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();

        assert!(loaded.is_transcript_keyed());
        assert_eq!(loaded.csq_columns(), ["pLI".to_string()]);
        assert_eq!(loaded.annotate_transcript("ENST00000318602.11").unwrap().unwrap(), r#"{"pLI":"1"}"#);
        assert_eq!(loaded.annotate_transcript("NM_130786.4").unwrap().unwrap(), r#"{"pLI":"2"}"#);
        assert!(loaded.annotate_transcript("ENST00000000001").unwrap().is_none());
    }

    #[test]
    fn symbol_keyed_file_is_not_transcript_keyed() {
        let header = GeneHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "omim".into(),
            name: "OMIM".into(),
            version: "1".into(),
            assembly: "GRCh38".into(),
        };
        let mut index = GeneIndex::new(header);
        index.add(GeneRecord { gene_symbol: "TP53".into(), json: "{}".into() });
        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        assert_eq!(&buf[..8], OGA_MAGIC);
        let loaded = GeneIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();
        assert!(!loaded.is_transcript_keyed());
        assert!(loaded.annotate_gene("TP53").unwrap().is_some());
    }
}
