# Greengenes2 ↔ MetaPhlAn genus lookup

The included crosswalk is ready to use with Python 3.9+. It uses only the
standard library. Keep `genus_crosswalk.tsv.gz` beside `gg2_metaphlan.py`.
Run from the extracted folder.

## Look up genera

Provide names as command arguments or put one genus per line in a text file.
Keep suffixes such as `_388775` and identifiers such as `GGB9278`.

```bash
python3 gg2_metaphlan.py lookup 16s-to-shotgun --input gg2_names.txt --output shotgun_names.tsv
python3 gg2_metaphlan.py lookup shotgun-to-16s --input shotgun_names.txt --output gg2_names.tsv
```

Try individual queries using the included crosswalk.

```bash
python3 gg2_metaphlan.py lookup 16s-to-shotgun Thomasclavelia --output example_sg.tsv
python3 gg2_metaphlan.py lookup shotgun-to-16s GGB9278 Erysipelatoclostridium --output example_16s.tsv
```

The `match` column uses the target database's original genus label.
Multiple candidate labels get separate rows. Support counts are distinct
SGB counts per pair. An SGB can support more than one pair, so the counts
are unsuitable as denominators for abundance or confidence estimates.
`no_supported_match` means the query genus is present in the source and the
available evidence produced zero links. `not_in_source_reference` means the
label did not occur in the source reference used for this build.

## How the bundled crosswalk was built

1. Read the species and genus labels from the supplied Greengenes2 2024.09
   V4 SINTAX FASTA.
2. Read genus and SGB membership from the taxonomy metadata of the original
   `mpa_vJun23_CHOCOPhlAnSGB_202307` index.
3. Link SGBs to GTDB R220 using the bioBakery Jun23 mapping. The mapping file
   uses a March 2024 release suffix. The June 2023 index update retained the
   SGB IDs and corrected some taxonomy labels. Joining by SGB ID preserves
   the original July 2023 genus labels on output.
4. Link a GTDB species to GG2 when the species label occurs in both or when
   the taxon-sync row lists a single species label on each side of a shared
   NCBI TaxID. The GG2 species must occur in the actual SINTAX FASTA.
   Genus pairs are derived from these species-level links.

These are reference-supported candidates. Matching species labels and
TaxIDs can reflect different taxon boundaries in the source databases.
The mapping includes no V4 sequence alignment or classifier rerun.
The compressed TSV records the evidence per SGB, the GTDB species and
genus, the GG2 species and any supporting NCBI TaxIDs.

## Rebuild

Use the original `.pkl` file from your trusted MetaPhlAn installation.
Pickle files can execute code when loaded. Use only trusted sources.
For a smaller build limited to observed taxa, replace `--metaphlan-db`
with `--shotgun-csv` and supply a CSV with MetaPhlAn lineage column headers.

```bash
curl -fL 'https://forum.biobakery.org/uploads/short-url/5RsfJ2AvwSpPp0wMWjxRUnCk00i.tsv' -o sgb2gtdb_r220.tsv
curl -fL 'https://www.taxonsync.org/api/v1/dataset?gtdb=R220&gg2=2024.09' -o taxonsync.tsv

python3 gg2_metaphlan.py build \
  --gg2-fasta gg2-2024.09-v4-sintax.fasta \
  --metaphlan-db mpa_vJun23_CHOCOPhlAnSGB_202307.pkl \
  --sgb2gtdb sgb2gtdb_r220.tsv \
  --taxonsync taxonsync.tsv \
  --output genus_crosswalk.tsv
gzip -nf genus_crosswalk.tsv
```

Sources

- [Greengenes2 2024.09 release](https://ftp.microbio.me/greengenes_release/2024.09/)
- [Original MetaPhlAn index archive](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.tar)
- [bioBakery Jun23 SGB-to-GTDB R220 assignment](https://forum.biobakery.org/t/update-for-sgb-to-gtdb-translation-file-r220/7388/4)
- [Merge compatibility of the Jun23 index releases](https://forum.biobakery.org/t/announcing-metaphlan-4-1-release-new-virome-database-and-sgb-database-update/6668)
- [taxon-sync export methodology](https://www.taxonsync.org/docs)
