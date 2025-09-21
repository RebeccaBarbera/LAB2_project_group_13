import sys
import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        # The regular expression is used to extract the next link for pagination
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

# This function actually retrieve the next data batch in the search.
# The function act as an iterator, yielding the next result batch at every call
# The function terminates after the last batch has been returned. In this case,
# the next link will be None
def get_batch(batch_url):
    while batch_url:
        # Run the API call
        response = session.get(batch_url)
        # Will raise an error if an error status code is obtained
        response.raise_for_status()
        # Get the total number of entries in the search
        total = response.headers["x-total-results"]
        # Yield the response and the total number of entries
        yield response, total
        # Get the link to the API call for the next data batch
        batch_url = get_next_link(response.headers)



def filter_entry(entry):
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

batch_size = 500
url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29+AND+%28length%3A%5B40+TO+*+%5D%29+AND+%28taxonomy_id%3A2759%29+NOT+%28ft_signal%3A*%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29+OR+%28cc_scl_term_exp%3ASL-0039%29%29%29&size=500"
# We set the name of the output file, we want TSV output
output_file = "filtered_negative.tsv"

# We define a function to better control the TSV format in output.
# In particular, we run the API call requiring JSON format and build our own TSV file
# The this aim, the following function extract and process specific fields from the JSON file

def extract_fields(entry):
    # We extract the accession, the sequence length and the start and end location of the signal peptide
    s, e = 0, 0
    # We iterate over the features of the entry
    for f in entry["features"]:
        # We only consider the SP segment
        if f["type"] == "Signal":
            e = f["location"]["end"]["value"]
            break
 
    organism_name = entry["organism"]["scientificName"]
    protein_length = entry["sequence"]["length"]
    accession = entry["primaryAccession"]
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
    return (accession, organism_name, kingdom, protein_length, e)


def get_dataset(search_url, filter_function, extract_function, output_file_name):
    filtered_json = []
    n_total, n_filtered = 0, 0
    # Run the API call in batches
    for batch, total in get_batch(search_url):
        # parse the JSON body of the response
        batch_json = json.loads(batch.text)
        # filter the entries
        for entry in batch_json["results"]:
            n_total += 1
            # Check if the entry passes the filter
            if filter_function(entry):
                n_filtered += 1
                filtered_json.append(extract_function(entry))
    print(n_total, n_filtered)
    with open(output_file_name, "w") as ofs:
        # Write header row
        print("UniprotAccession\tOrganism\tKingdom\tProteinLength\tSPPosition", file=ofs)
        for entry in filtered_json:
            # Print the fields in TSV format
            print(*entry, sep="\t", file=ofs)
        ofs.close

get_dataset(url, filter_entry, extract_fields, output_file)
