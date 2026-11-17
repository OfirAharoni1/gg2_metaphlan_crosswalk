#!/usr/bin/env python3
"""Build and query a Greengenes2 V4 / MetaPhlAn Jun23 genus crosswalk."""

import argparse
import bz2
import csv
import gzip
import pickle
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


COLUMNS = ("gg2_genus", "metaphlan_genus", "sgb", "gtdb_genus", "gtdb_species",
           "gg2_species", "evidence", "ncbi_taxids", "status")
SGB = re.compile(r"SGB\d+(?:_group)?")
GENUS_PREFIX = re.compile(r"^g(?::|__)", re.IGNORECASE)


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def clean_label(value):
    return " ".join(value.split()).casefold()


def parse_lineage(lineage, separator, rank_separator):
    ranks = {}
    for token in lineage.split(separator):
        rank, found, label = token.strip().partition(rank_separator)
        if found:
            ranks[rank] = label.strip()
    return ranks


def gg2_reference(path):
    genera = set()
    species = defaultdict(set)
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith(">") or ";tax=" not in line:
                continue
            ranks = parse_lineage(line.split(";tax=", 1)[1].rstrip().rstrip(";"), ",", ":")
            genus = ranks.get("g")
            if genus and ranks.get("d") in {"Bacteria", "Archaea"}:
                genera.add(genus)
                if ranks.get("s"):
                    species[clean_label(ranks["s"])].add((genus, ranks["s"]))
    if not genera:
        raise ValueError("FASTA contains no SINTAX bacterial or archaeal genus assignments")
    return genera, species


def gtdb_reference(path):
    assignments = defaultdict(set)
    with open_text(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 2 or not SGB.fullmatch(row[0]):
                raise ValueError(f"Invalid SGB2GTDB row: {row!r}")
            ranks = parse_lineage(row[1], ";", "__")
            assignments[row[0]].add((ranks.get("g", ""), ranks.get("s", "")))
    if not assignments:
        raise ValueError("SGB2GTDB table is empty")
    return assignments


def taxonsync_links(path):
    links = defaultdict(lambda: defaultdict(set))
    with open_text(path) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        expected = {"ncbi_taxid", "gtdb_R220", "gg2_2024.09"}
        if not expected.issubset(rows.fieldnames or []):
            raise ValueError("Expected a taxon-sync dataset pinned to GTDB R220 and GG2 2024.09")
        for row in rows:
            gtdb, gg2 = row["gtdb_R220"], row["gg2_2024.09"]
            if gtdb and gg2 and " / " not in gtdb and " / " not in gg2:
                links[clean_label(gtdb)][clean_label(gg2)].add(row["ncbi_taxid"])
    return links


def shotgun_lineages(args):
    if args.metaphlan_db:
        # MetaPhlAn stores its trusted reference metadata as a bzip2 pickle.
        with bz2.open(args.metaphlan_db, "rb") as handle:
            return list(pickle.load(handle)["taxonomy"])
    if args.shotgun_csv:
        with open_text(args.shotgun_csv) as handle:
            return next(csv.reader(handle))
    with open_text(args.shotgun_lineages) as handle:
        return [line.strip() for line in handle if line.strip()]


def shotgun_reference(args):
    genera = set()
    members = set()
    for lineage in shotgun_lineages(args):
        separator = "|" if "|" in lineage else ","
        ranks = parse_lineage(lineage, separator, "__")
        if ranks.get("k", ranks.get("d")) not in {"Bacteria", "Archaea"}:
            continue
        genus = ranks.get("g")
        if genus:
            genera.add(genus)
            sgb = ranks.get("t", "")
            if SGB.fullmatch(sgb):
                members.add((genus, sgb))
    if not members:
        raise ValueError("Shotgun source contains no bacterial or archaeal g__/t__SGB lineages")
    return genera, members


def build_rows(gg_genera, gg_species, mpa_genera, members, assignments, links):
    for genus in sorted(gg_genera):
        yield genus, "", "", "", "", "", "", "", "present_gg2_reference"
    for genus in sorted(mpa_genera):
        yield "", genus, "", "", "", "", "", "", "present_metaphlan_source"
    for mpa_genus, sgb in sorted(members):
        if sgb not in assignments:
            yield "", mpa_genus, sgb, "", "", "", "", "", "no_r220_assignment"
            continue
        for gtdb_genus, gtdb_species in sorted(assignments[sgb]):
            key = clean_label(gtdb_species)
            evidence = defaultdict(set)
            if key:
                for pair in gg_species.get(key, ()):
                    evidence[pair]  # Exact species-label match.
                for gg_name, ids in links.get(key, {}).items():
                    for pair in gg_species.get(gg_name, ()):
                        evidence[pair].update(ids)
            if not evidence:
                status = "gtdb_genus_only" if gtdb_genus and not key else "no_gg2_species_link"
                if not key and not gtdb_genus:
                    status = "no_r220_species_or_genus"
                yield "", mpa_genus, sgb, gtdb_genus, gtdb_species, "", "", "", status
                continue
            for (gg_genus, gg_name), taxids in sorted(evidence.items()):
                methods = []
                if key == clean_label(gg_name):
                    methods.append("shared_species_label")
                if taxids:
                    methods.append("taxonsync_single_name_taxid")
                yield (gg_genus, mpa_genus, sgb, gtdb_genus, gtdb_species, gg_name,
                       ",".join(methods), ",".join(sorted(taxids)), "candidate")


def build(args):
    gg_genera, gg_species = gg2_reference(args.gg2_fasta)
    mpa_genera, members = shotgun_reference(args)
    assignments = gtdb_reference(args.sgb2gtdb)
    links = taxonsync_links(args.taxonsync) if args.taxonsync else {}
    counts = Counter()
    pairs = set()
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in build_rows(gg_genera, gg_species, mpa_genera, members, assignments, links):
            writer.writerow(row)
            counts[row[-1]] += 1
            if row[-1] == "candidate":
                pairs.add(row[:2])
    if not pairs:
        raise ValueError("No candidate pairs. Check input formats and reference releases.")
    print(f"{len(members)} genus/SGB memberships, {len(pairs)} genus pairs; {dict(counts)}",
          file=sys.stderr)


def lookup(args):
    source = "gg2_genus" if args.direction == "16s-to-shotgun" else "metaphlan_genus"
    target = "metaphlan_genus" if source == "gg2_genus" else "gg2_genus"
    index = defaultdict(list)
    with open_text(args.map) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not set(COLUMNS).issubset(rows.fieldnames or []):
            raise ValueError("Invalid crosswalk header")
        for row in rows:
            if row[source]:
                index[clean_label(row[source])].append(row)
    queries = list(args.names)
    if args.input:
        with open_text(args.input) as handle:
            queries.extend(line.strip() for line in handle if line.strip())
    if not queries:
        raise ValueError("Enter genus names or supply --input")
    fields = ("query", "match", "status", "supporting_sgbs", "supporting_sgb_count",
              "gtdb_genera", "evidence", "unlinked_source_sgbs")
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for query in dict.fromkeys(queries):
            rows = index.get(clean_label(GENUS_PREFIX.sub("", query.strip())))
            matches = defaultdict(list)
            for row in rows or ():
                if row["status"] == "candidate":
                    matches[row[target]].append(row)
            unlinked = set()
            if source == "metaphlan_genus":
                unlinked = {row["sgb"] for row in rows or () if row["sgb"]} - {
                    row["sgb"] for match_rows in matches.values() for row in match_rows
                }
            if not matches:
                status = "no_supported_match" if rows else "not_in_source_reference"
                writer.writerow((query, "", status, "", 0, "", "", ",".join(sorted(unlinked))))
            for label, evidence in sorted(matches.items()):
                sgbs = {row["sgb"] for row in evidence}
                gtdb = {row["gtdb_genus"] for row in evidence if row["gtdb_genus"]}
                methods = {method for row in evidence for method in row["evidence"].split(",")}
                status = "multiple_candidates" if len(matches) > 1 else "candidate"
                writer.writerow((query, label, status, ",".join(sorted(sgbs)), len(sgbs),
                                 ",".join(sorted(gtdb)), ",".join(sorted(methods)),
                                 ",".join(sorted(unlinked))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    builder = sub.add_parser("build", help="Build a version-specific crosswalk once")
    builder.add_argument("--gg2-fasta", required=True, type=Path)
    builder.add_argument("--sgb2gtdb", required=True, type=Path)
    builder.add_argument("--taxonsync", type=Path)
    reference = builder.add_mutually_exclusive_group(required=True)
    reference.add_argument("--metaphlan-db", type=Path, help="Trusted MetaPhlAn .pkl")
    reference.add_argument("--shotgun-csv", type=Path, help="Formatted table with taxonomy headers")
    reference.add_argument("--shotgun-lineages", type=Path, help="One MetaPhlAn lineage per line")
    builder.add_argument("--output", required=True, type=Path)
    finder = sub.add_parser("lookup", help="Look up genus labels in either direction")
    finder.add_argument("direction", choices=("16s-to-shotgun", "shotgun-to-16s"))
    finder.add_argument("names", nargs="*")
    finder.add_argument("--input", type=Path, help="One query genus per line")
    finder.add_argument("--map", type=Path,
                        default=Path(__file__).resolve().with_name("genus_crosswalk.tsv.gz"))
    finder.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    inputs = [value for key, value in vars(args).items()
              if key != "output" and isinstance(value, Path)]
    if args.output.resolve() in {path.resolve() for path in inputs}:
        parser.error("output must differ from every input")
    try:
        (build if args.command == "build" else lookup)(args)
    except (OSError, ValueError, csv.Error, KeyError, pickle.UnpicklingError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
