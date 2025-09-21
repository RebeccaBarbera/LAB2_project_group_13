import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re

# sessione con retry
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500,502,503,504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# funzione per paginazione
def get_next_link(headers):
    if "Link" in headers:
        match = re.match(r'<(.+)>; rel="next"', headers["Link"])
        if match:
            return match.group(1)

# iterator per batch
def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        yield response, response.headers.get("x-total-results", "unknown")
        batch_url = get_next_link(response.headers)

# estrazione campi positivo
def extract_fields(entry):
    accession = entry["primaryAccession"]
    organism_name = entry["organism"]["scientificName"]
    protein_length = entry["sequence"]["length"]
    kingdom = "Other"
    for item in entry["organism"]["lineage"]:
        if item == "Metazoa":
            kingdom = "Metazoa"
            break
        elif item == "Fungi":
            kingdom = "Fungi"
            break
        elif item == "Viridiplantae":
            kingdom = "Viridiplantae"
            break

    return accession, organism_name, kingdom, protein_length

# filtro positivo
def filter_positive(entry):
     # We iterate over the features of the entry
    for feature in entry["features"]:
        # We only consider features of type Signal
        if feature["type"] == "Signal":
            # Check if the signal peptide is longer tha 14 residues
            end = feature["location"]["end"]["value"]
            start = feature["location"]["start"]["value"]
            cleaved =  feature["description"]
                
            if start is None or end is None:
                    continue

            if end >= 14 and cleaved == "":

                return True
    return False

# funzione negativa + estrazione
def neg_entries(entry):
    accession, organism_name, kingdom, protein_length = extract_fields(entry)
    tm_first90 = False
    for feature in entry["features"]:
        if feature["type"] == "Transmembrane":
            start = feature["location"]["start"]["value"]
            if start <= 90:
                tm_first90 = True
                break
    return accession, organism_name, kingdom, protein_length, tm_first90

# funzione per creare TSV
def get_positive(search_url, filter_function, extract_function, output_file_name, header):
    n_total, n_filtered = 0,0
    with open(output_file_name, "w") as ofs:
        print(header, file=ofs)
        for batch in get_batch(search_url):
            response, _ = batch  # Unpack tuple
            batch_json = json.loads(response.text)
            for entry in batch_json.get("results", []):
                n_total += 1
                if filter_function(entry):
                    n_filtered += 1
                    fields = extract_function(entry)
                    print(*fields, sep="\t", file=ofs)
    print(f"Total entries: {n_total}, Filtered entries: {n_filtered}, \n writing results in {output_file_name}...")

def get_negative(search_url, extract_function, output_file_name, header):
    n_total = 0
    with open(output_file_name, "w") as ofs:
        print(header, file=ofs)
        for batch in get_batch(search_url):
            response, _ = batch  # Unpack tuple
            batch_json = json.loads(response.text)
            for entry in batch_json.get("results", []):
                n_total += 1
                fields = extract_function(entry)
                print(*fields, sep="\t", file=ofs)
    print(f"Total entries: {n_total}, \n writing results in {output_file_name}....")

# URL positivo e negativo
positive_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)%29&size=500"
negative_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (reviewed:true) AND (existence:1) AND (length:[40 TO *]) AND (taxonomy_id:2759) NOT (ft_signal:*) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))%29&size=500"

# dataset positivo
get_positive(positive_url, filter_positive, extract_fields,
            "positive_set.tsv",
            "UniprotAccession\tOrganism\tKingdom\tProteinLength\tSPPosition")

# dataset negativo
get_negative(negative_url, neg_entries,
            "negative_set.tsv",
            "UniprotAccession\tOrganism\tKingdom\tProteinLength\tTMHelixFirst90")
