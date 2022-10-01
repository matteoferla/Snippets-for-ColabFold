# Pyrosetta code will be integrated in future into other module

from Bio import pairwise2
from IPython.display import HTML
import requests
import pyrosetta_help as ph
import functools


def retrieve_uniprot_data(uniprot: str) -> dict:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.json')
    response.raise_for_status()
    return response.json()

# moved over to PH.
fix_offset = ph.fix_offset
assign_chain_letter = ph.assign_chain_letter
what_is_chain = ph.what_is_chain
