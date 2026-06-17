#!/usr/bin/env python3
# fetch_epsps_references.py - fetch one representative EPSPS (aroA) protein
# per genus from NCBI, building refs/epsps/reference_epsps.faa automatically.
#
# Input: results/r/epsps_genera_to_curate.csv (written by 05_epsps_overlay.R)
# Output: refs/epsps/reference_epsps.faa, header format >Genus|accession
#         results/epsps/fetch_log.csv, one row per genus with the outcome
#
# Requires NCBI_API_KEY as an environment variable. Without it the script
# still runs, at the unkeyed rate limit of 3 requests/second instead of 10.
#
# Search strategy: for each genus, query the NCBI protein database for
# aroA or EPSP synthase annotated entries restricted to that genus, prefer
# a RefSeq hit if one exists, otherwise take the first hit. Genera with no
# hit are logged as not_found, not silently dropped.

import csv
import os
import re
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

PROJECT = os.path.expanduser("~/pilot2018")
GENERA_CSV = f"{PROJECT}/results/r/epsps_genera_to_curate.csv"
OUT_FAA = f"{PROJECT}/refs/epsps/reference_epsps.faa"
LOG_CSV = f"{PROJECT}/results/epsps/fetch_log.csv"

API_KEY = os.environ.get("NCBI_API_KEY", "")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DELAY = 0.11 if API_KEY else 0.34

def looks_like_silva_clade_code(name):
    if re.search(r'(clade|lineage|marine_group)$', name, re.IGNORECASE):
        return True
    if re.match(r'^[A-Za-z]{1,4}[\d]', name) and any(c.isdigit() for c in name):
        return True
    if re.match(r'^[A-Z0-9]+[\-_][A-Z0-9]+$', name) and not re.match(r'^[A-Za-z]+-[A-Za-z]+$', name):
        return True
    return False

def eutils_get(endpoint, params):
    params = dict(params)
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()

def search_genus(genus):
    name = genus.replace("g__", "").strip()
    if not name or name.lower() in ("unknown", "unknown_family", "uncultured", "na", "metagenome"):
        return None, "skipped_placeholder", name
    if looks_like_silva_clade_code(name):
        return None, "skipped_silva_clade_code", name

    candidates_names = [name]
    if "-" in name and not name.startswith("-"):
        parts = name.split("-")
        if len(parts) == 2 and all(p.isalpha() for p in parts):
            candidates_names = parts + [name]

    for query_name in candidates_names:
        query = f'(aroA[Gene Name] OR "EPSP synthase"[Title]) AND "{query_name}"[Organism]'
        xml = eutils_get("esearch.fcgi", {
            "db": "protein", "term": query, "retmax": 20, "sort": "relevance"
        })
        root = ET.fromstring(xml)
        ids = [e.text for e in root.findall(".//Id")]
        if ids:
            return ids, "found", query_name
        time.sleep(DELAY)

    return None, "not_found", name

def fetch_summaries(ids):
    xml = eutils_get("esummary.fcgi", {
        "db": "protein", "id": ",".join(ids), "retmode": "xml"
    })
    root = ET.fromstring(xml)
    summaries = []
    for doc in root.findall(".//DocSum"):
        uid = doc.findtext("Id")
        caption = None
        for item in doc.findall("Item"):
            if item.get("Name") == "Caption":
                caption = item.text
        summaries.append((uid, caption or ""))
    return summaries

def pick_best(ids):
    summaries = fetch_summaries(ids)
    refseq = [s for s in summaries if s[1].startswith(("WP_", "NP_", "YP_"))]
    chosen = refseq[0] if refseq else summaries[0]
    return chosen

def fetch_fasta(uid):
    return eutils_get("efetch.fcgi", {
        "db": "protein", "id": uid, "rettype": "fasta", "retmode": "text"
    })

def load_previous_results():
    prev_faa_records = []
    prev_ok_genera = set()
    if os.path.exists(LOG_CSV) and os.path.exists(OUT_FAA):
        with open(LOG_CSV) as f:
            for row in csv.DictReader(f):
                if row["status"] == "ok":
                    prev_ok_genera.add(row["genus"])
        with open(OUT_FAA) as f:
            prev_faa_records = [l for l in f]
    return prev_faa_records, prev_ok_genera

def main():
    if not os.path.exists(GENERA_CSV):
        sys.exit(f"Genera list not found at {GENERA_CSV}. Run 05_epsps_overlay.R first.")

    os.makedirs(os.path.dirname(OUT_FAA), exist_ok=True)
    os.makedirs(os.path.dirname(LOG_CSV), exist_ok=True)

    with open(GENERA_CSV) as f:
        reader = csv.DictReader(f)
        genera = [row["Genus"] for row in reader]

    prev_faa_records, prev_ok_genera = load_previous_results()
    total_genera = len(genera)
    to_fetch = [g for g in genera if g not in prev_ok_genera]

    print(f"{total_genera} genera total. {len(prev_ok_genera)} already fetched in a prior run.")
    print(f"Fetching the remaining {len(to_fetch)} genera. This retries any not_found")
    print("entries from a prior run with improved compound-name and clade-code handling.")
    print(f"NCBI API key {'set' if API_KEY else 'not set, running at the unkeyed rate limit'}.")

    faa_records = list(prev_faa_records)
    log_rows = []
    if os.path.exists(LOG_CSV):
        with open(LOG_CSV) as f:
            log_rows = [row for row in csv.DictReader(f) if row["status"] == "ok"]

    genera_to_process = to_fetch

    for i, genus in enumerate(genera_to_process, 1):
        try:
            ids, status, matched_name = search_genus(genus)
            time.sleep(DELAY)
            if status != "found":
                log_rows.append({"genus": genus, "status": status, "accession": "", "matched_name": ""})
                print(f"[{i}/{len(genera_to_process)}] {genus}: {status}")
                continue

            uid, accession = pick_best(ids)
            time.sleep(DELAY)
            fasta = fetch_fasta(uid)
            time.sleep(DELAY)

            if not fasta.strip().startswith(">"):
                log_rows.append({"genus": genus, "status": "fetch_failed", "accession": accession, "matched_name": matched_name})
                print(f"[{i}/{len(genera_to_process)}] {genus}: fetch_failed")
                continue

            seq_lines = fasta.strip().splitlines()
            seq = "".join(seq_lines[1:])
            clean_genus = genus.replace("g__", "").strip().replace(" ", "_")
            faa_records.append(f">{clean_genus}|{accession}\n{seq}\n")
            log_rows.append({"genus": genus, "status": "ok", "accession": accession, "matched_name": matched_name})
            note = f" (matched as {matched_name})" if matched_name != genus.replace("g__", "").strip() else ""
            print(f"[{i}/{len(genera_to_process)}] {genus}: ok ({accession}){note}")

        except Exception as e:
            log_rows.append({"genus": genus, "status": f"error: {e}", "accession": "", "matched_name": ""})
            print(f"[{i}/{len(genera_to_process)}] {genus}: error: {e}")
            time.sleep(DELAY)

    with open(OUT_FAA, "w") as f:
        f.writelines(faa_records)

    with open(LOG_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["genus", "status", "accession", "matched_name"])
        w.writeheader()
        w.writerows(log_rows)

    ok = sum(1 for r in log_rows if r["status"] == "ok")
    print(f"\nDone. {ok}/{total_genera} genera fetched successfully overall.")
    not_found_n = sum(1 for r in log_rows if r["status"] == "not_found")
    clade_n = sum(1 for r in log_rows if r["status"] == "skipped_silva_clade_code")
    print(f"{not_found_n} not_found (real genus name, no NCBI hit).")
    print(f"{clade_n} skipped_silva_clade_code (SILVA placeholder clade label, not a queryable organism).")
    print(f"Reference FASTA: {OUT_FAA}")
    print(f"Fetch log (review failures here): {LOG_CSV}")

if __name__ == "__main__":
    main()
