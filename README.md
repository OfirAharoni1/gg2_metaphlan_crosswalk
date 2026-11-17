# Greengenes2 ↔ MetaPhlAn genus crosswalk

Look up candidate matches between genus labels from 16S (Greengenes2)
and shotgun metagenomics (MetaPhlAn), in either direction.

Requires Python 3.9+. No packages to install.

## Quick start

Create `query_16s.txt` with one genus label per line. Create
`shotgun_reference.csv` with a `genus` header and one target genus per row.

```bash
python3 genus_crosswalk.py --s16 query_16s.txt \
  --shotgun shotgun_reference.csv --output matches.tsv
```

## Example output

This example uses `Bacteroides`, `ExampleGenus`, and `uncultured` as 16S
queries. The target reference contains `Bacteroides`.

| label_16s | label_sg | status | evidence | ncbi_genus_ids_16s | ncbi_genus_ids_sg | genus_sg | genus_16s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bacteroides | Bacteroides | `exact_label_candidate` | | | | Bacteroides | Bacteroides |
| ExampleGenus | | `unmatched_16s` | | | | | ExampleGenus |
| uncultured | | `unresolved_16s` | | | | | |

### Output columns

| Column | Meaning |
| --- | --- |
| `label_16s` | Original label or taxonomic lineage on the 16S side. |
| `label_sg` | Original label or taxonomic lineage on the shotgun side. |
| `status` | Candidate match, curated relationship, or reason no match was reported. |
| `evidence` | Source or justification supplied for a curated pair. |
| `ncbi_genus_ids_16s` | NCBI genus TaxIDs associated with the 16S name when `--taxdump` is used. Multiple IDs are separated by semicolons. |
| `ncbi_genus_ids_sg` | NCBI genus TaxIDs associated with the shotgun name when `--taxdump` is used. Multiple IDs are separated by semicolons. |
| `genus_sg` | Genus name extracted from `label_sg`. |
| `genus_16s` | Genus name extracted from `label_16s`. |

An exact-label match means the genus names agree after case and whitespace
normalization. Treat automatic matches as candidates for taxonomic equivalence.
A query can produce several rows, one per candidate. `unmatched` means the
tool extracted a genus but found no candidate in the supplied reference.
`unresolved` means it could not extract a usable genus name.

## Shotgun → 16S

```bash
python3 genus_crosswalk.py --sg-query query_sg.txt \
  --16s-reference 16s_reference.csv --output matches.tsv
```

The output columns keep the same order in both directions.

<details>
<summary>Other inputs and matching options</summary>

- **Query tables.** CSV and TSV files with a `genus` column also work.
  Use `--s16-column` or `--shotgun-column` for a different column name.
- **Shotgun abundance table.** Pass `--shotgun-profile profile.csv` in
  place of `--shotgun` to read genus-level taxonomy from column headers.
- **Default references.** If you omit the reference flag, the script
  looks for `shotgun_reference.csv` or `16s_reference.csv` next to the
  script, depending on the direction.
- **Reference coverage.** A reference built from your samples can only
  return genera present in those samples. For database-wide lookup,
  use genus labels from the release used by your classifier.
- **NCBI names and synonyms.** Pass `--taxdump /path/to/taxdump` with
  `nodes.dmp` and `names.dmp`. The script can report
  `shared_ncbi_genus_candidate` or flag ambiguity as `ambiguous_candidate`.
- **Curated pairs.** Pass `--reviewed reviewed.tsv` with columns
  `genus_16s`, `genus_sg`, `relation`, and `evidence`. The relation
  must be `equivalent`, `partial`, or `distinct`; it appears in the output
  as `reviewed_equivalent`, `reviewed_partial`, or `reviewed_distinct`.

</details>