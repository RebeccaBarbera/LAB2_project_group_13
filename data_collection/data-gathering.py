import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re

# Setup requests session with retry logic
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500,502,503,504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re.match(r'<(.+)>; rel="next"', headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        yield response, response.headers.get("x-total-results", "unknown")
        batch_url = get_next_link(response.headers)

def extract_fields(entry, sp_position=None):
    organism_name = re.sub(r'\s*\(.*\)$', '', entry["organism"]["scientificName"])
    accession = entry["primaryAccession"]
    protein_length = entry["sequence"]["length"]
    kingdom = "Other"
    for item in entry["organism"]["lineage"]:
        if item == "Metazoa":
            kingdom = "Metazoa"
        elif item == "Fungi":
            kingdom = "Fungi"
        elif item == "Viridiplantae":
            kingdom = "Viridiplantae"
    return accession, organism_name, kingdom, protein_length, sp_position

def filter_positive(entry):
    # Returns True and SPPosition if a valid signal peptide is found
    for feature in entry["features"]:
        if feature["type"] == "Signal":
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            cleaved = feature["description"]
            if start is None or end is None:
                continue
            if end >= 14 and cleaved == "":
                return True, f"{end}"
    return False, None

def neg_entries(entry):
    accession, organism_name, kingdom, protein_length, _ = extract_fields(entry, sp_position=None)
    tm_first90 = False
    for feature in entry["features"]:
        if feature["type"] == "Transmembrane":
            start = feature["location"]["start"]["value"]
            if start <= 90:
                tm_first90 = True
                break
    return accession, organism_name, kingdom, protein_length, tm_first90

def get_positive(search_url, filter_function, extract_function, output_tsv, output_fasta, header):
    n_total, n_filtered = 0, 0
    with open(output_tsv, "w") as ofs, open(output_fasta, "w") as fasta_ofs:
        print(header, file=ofs)
        for batch in get_batch(search_url):
            response, _ = batch
            batch_json = json.loads(response.text)
            for entry in batch_json.get("results", []):
                n_total += 1
                passed, sp_position = filter_function(entry)
                if passed:
                    n_filtered += 1
                    fields = extract_function(entry, sp_position)
                    print(*fields, sep="\t", file=ofs)
                    seq = entry["sequence"]["value"]
                    fasta_ofs.write(f">{fields[0]}\n{seq}\n")
    print(f"Total entries: {n_total}, Filtered entries: {n_filtered}, \nwriting results in {output_tsv} and {output_fasta}...")

def get_negative(search_url, extract_function, output_tsv, output_fasta, header):
    n_total = 0
    with open(output_tsv, "w") as ofs, open(output_fasta, "w") as fasta_ofs:
        print(header, file=ofs)
        for batch in get_batch(search_url):
            response, _ = batch
            batch_json = json.loads(response.text)
            for entry in batch_json.get("results", []):
                n_total += 1
                fields = extract_function(entry)
                print(*fields, sep="\t", file=ofs)
                seq = entry["sequence"]["value"]
                fasta_ofs.write(f">{fields[0]}\n{seq}\n")
    print(f"Total entries: {n_total}, \nwriting results in {output_tsv} and {output_fasta}...")

# UniProt API queries for positive and negative datasets
positive_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)%29&size=500"
negative_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (reviewed:true) AND (existence:1) AND (length:[40 TO *]) AND (taxonomy_id:2759) NOT (ft_signal:*) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))%29&size=500"

get_positive(
    positive_url, filter_positive, extract_fields,
    "positive_set.tsv",
    "positive_set.fasta",
    "UniprotAccession\tOrganism\tKingdom\tProteinLength\tSPPosition"
)

get_negative(
    negative_url, neg_entries,
    "negative_set.tsv",
    "negative_set.fasta",
    "UniprotAccession\tOrganism\tKingdom\tProteinLength\tTMHelixFirst90"
)