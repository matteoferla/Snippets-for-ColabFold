"""
Chop up an AlphaFold PDB file into multiple PDB files, based on contiguous chains with pDTTL >= cutoff
"""

import re
import os
import requests


def split_n_save(uniprot: str, folder: str = '.') -> i:
    """
    Given a uniprot accession, chop it up and write it to a folder
    :param uniprot: the gibberish id (accession), not the name-like id (Uniprot name)
    :param folder:
    :return: number of files
    """
    if not os.path.exists(folder):
        os.mkdir(folder)
    master_pdb_block = get_pdbblock(uniprot)
    pdb_blocks = split_by_LDDT(master_pdb_block)
    for i, pdb_block in enumerate(pdb_blocks):
        with open(os.path.join(folder, f'{uniprot}.{i}.pdb'), 'w') as fh:
            fh.write(pdb_block)
    return i + 1


def get_pdbblock(uniprot: str, version=3) -> str:
    """
    Retrieve the PDB block for a given uniprot accession.
    Note that the version ticks up once a year or so.
    1 = release of key organism SwissProt entries
    2 = release of all SwissProt entries
    3 = release of all Trembl entries

    :param uniprot:
    :param version:
    :return:
    """
    response: requests.Response = requests.get(
        f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v{version}.pdb')
    response.raise_for_status()
    return response.text


def split_by_LDDT(pdbblock: str, cutoff: float = 70):
    """
    Split barbarically the block by the LDDT score.
    Will not work for HETATM and partial occupancy etc etc.
    It's for AlphaFold as of 2022.

    :param pdbblock:
    :param cutoff:
    :return:
    """
    subblock = ''
    subblocks = []
    serial = 0
    for line in pdbblock.split('\n'):
        if line[:4] != 'ATOM':
            continue
        if float(line[60:67].strip()) >= cutoff:
            serial += 1
            subblock += line[:6] + f'{serial: >5}' + line[11:] + '\n'
        elif subblock:
            subblocks.append(subblock)
            subblock = ''
            serial = 0
        else:
            pass
    if subblock:
        subblocks.append(subblock)
    return subblocks
