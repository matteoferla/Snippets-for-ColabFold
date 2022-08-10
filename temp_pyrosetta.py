# Pyrosetta code will be integrated in future into other module

from Bio import pairwise2
from IPython.display import display, HTML
import requests
import pyrosetta_help as ph
import functools


def retrieve_uniprot_data(uniprot: str) -> dict:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.json')
    response.raise_for_status()
    return response.json()


def fix_offset(pose, chain_i: int, ref_sequence: str) -> HTML:
    """
    Fix the pose numbering on ``chain_i`` based on ``ref_sequence``
    No gaps accepted in ref_sequence, so pad with Xs appropriately.
    """
    pdb_info = pose.pdb_info()
    assert not pdb_info.obsolete(), "pdb info is marked as obsolete"
    # ## Align
    alignment = pairwise2.align.globalms(ref_sequence,
                                         pose.chain_sequence(chain_i),
                                         1, -2,  # penalise mismatches by x2
                                         -1, -0.5  # penalise extension by 1/2
                                         )[0]
    assert alignment.seqA.find('-') == -1
    resi_map = [i + 1 for i, l in enumerate(alignment.seqB) if l != '-']
    pdb_info = pose.pdb_info()
    chain_resi_offset = -1
    for resi in ph.pose_range(pose):
        residue = pose.residue(resi)
        if residue.chain() != chain_i:
            continue
        if chain_resi_offset == -1:
            chain_resi_offset = resi - 1
        pdb_info.number(resi, resi_map[resi - chain_resi_offset - 1])
    # ## Output
    formatted = pairwise2.format_alignment(*alignment)
    a, gap, b, score = formatted.strip().split('\n')
    gap = ''.join(['.' if c == '|' else '*' for c in gap])
    return HTML(f'<div style="font-family:monospace; display: inline-block; white-space: nowrap;">' +
                f'{a}<br/>{gap.replace(" ", "*")}<br/>{b}<br/>{score}</div>')


def assign_chain_letter(pose, chain_i: int, new_letter: str):
    """
    Assign a new letter to chain number ``chain_i``
    """
    pdb_info = pose.pdb_info()
    assert not pdb_info.obsolete(), "pdb info is marked as obsolete"
    for resi in ph.pose_range(pose):
        residue = pose.residue(resi)
        if residue.chain() != chain_i:
            continue
        pdb_info.chain(resi, new_letter)



@functools.singledispatch
def what_is_chain(pose, chain_i: int):
    for resi in ph.pose_range(pose):
        residue = pose.residue(resi)
        if residue.chain() != chain_i:
            continue
        pdb_info = pose.pdb_info()
        return pdb_info.chain(resi)
    raise ValueError(f'No chain {chain_i}')


@what_is_chain.register
def _(pose, chain_name: str):
    pdb_info = pose.pdb_info()
    for resi in ph.pose_range(pose):
        if pdb_info.chain(resi) != chain_name:
            continue
        residue = pose.residue(resi)
        return residue.chain()
    raise ValueError(f'No chain {chain_name}')