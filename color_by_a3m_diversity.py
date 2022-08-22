from typing import List, Set, Callable

import functools
import operator
import os
import pymol2
import re
from Bio import SeqIO
from IPython.display import display, FileLink


def get_diversity(alignment_filename: str, gene_i: int = 0) -> List[Set[str]]:
    '''
    Returns a list of sets of AAs (with no gaps)
    '''
    # ## Parse
    assert os.path.splitext(alignment_filename)[1] == '.a3m', 'Only A3M files are supported'
    headers: List[str] = []
    seqs: List[str] = []
    with open(alignment_filename, 'r') as fh:
        _first_line = next(fh)[1:].strip().split()
        ranges: List[int] = list(map(int, _first_line[0].split(',')))
        cardinality: List[int] = list(map(int, _first_line[1].split(',')))
        for record in SeqIO.parse(fh, "fasta"):
            headers.append(record.description.split())
            seqs.append(str(record.seq).replace('\n', ''))
    # ## Process
    # the ranges need to be cumulative
    p = [0] + [sum(ranges[:i + 1]) for i in range(len(ranges))]
    cumranges = list(zip(p[:-1], p[1:]))
    # A3M has weird lowercases to mark insertions
    remove_weirdos: Callable[[str, ], str] = functools.partial(re.sub, '[a-z]', '')
    split_seqs: List[List[str]] = [[remove_weirdos(seq)[slice(*s)] for s in cumranges] for seq in seqs]
    gene_seqs: List[str] = list(map(operator.itemgetter(gene_i), split_seqs))
    # list of non-gap AAs per position
    return [set(map(operator.itemgetter(i), gene_seqs)).difference('-') for i in range(len(gene_seqs[0]))]


def map_diversity(pdbblock: str, diversity: List[Set[str]], outfile='colored.pse') -> None:
    """
    Colour the diversity from a A3M file as b factors in PyMol
    ```python
    pdbblock:str = requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v3.pdb').text
    diversity:List[Set[str]] = get_diversity('foo.a3m')
    map_diversity(pdbblock, diversity, outfile='foo.pse')
    ```
    """
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(pdbblock, 'alphafold')
        for i, d in enumerate(diversity):
            pymol.cmd.alter(f'alphafold and resi {i + 1}', f'b={len(d)}')
        pymol.cmd.sort()
        pymol.cmd.spectrum('b', 'white yellow red', minimum=1, maximum=10)
        pymol.cmd.save(outfile)
    display(FileLink(outfile))
