from __future__ import annotations
import functools
import io
import operator
import re
import requests
from typing import Set, List

import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import Record


class Novelty:

    def __init__(self):
        self.seen = set()

    def __call__(self, values):
        novel = len(self.seen.intersection(set(values))) == 0
        self.seen.update(values)
        return novel


def clean(description: str, **alternatives) -> str:
    """
    Strip fluff of the protein names.
    Additionally, passing a key value pair will replace the key with the value
    (remember ``foo(**{'hello world': 'bye Mars'})`` is okay, but ``foo(hello world='bye Mars')`` is not)
    """
    if description is None:
        description = ''
    description = str(description)
    for query in (r'LOW QUALITY PROTEIN:\s+', r' isoform \w+', r'PREDICTED:\s+',
                  'putative',
                  'uncharacterized protein', ', partial', '-like'):
        alternatives[query] = '' if query not in alternatives else str(alternatives[query])
    alternatives['  '] = ' '
    alternatives['LOC\d+'] = 'unknown'
    for query in alternatives:
        description = re.sub(query, alternatives[query], description)
    return description.strip()


def get_common_species(*dfs) -> Set[str]:
    """
    The dataframes (``Blaster.results``) contain a column called
    ``score`` (int) and ``specieses`` (the latter is a list of species names)
    """
    dfs: List[pd.DataFrame] = [df.loc[df.score > 50] for df in dfs]
    # set.union & set.intersection do return set
    species_sets: List[Set[str]] = [functools.reduce(set.union, df.specieses, set()) for df in dfs]  # noqa
    return functools.reduce(set.intersection, species_sets[1:], species_sets[0]) # noqa


class Blaster:

    def __init__(self, sequence: str):
        blast: io.StringIO = NCBIWWW.qblast(program='blastp',
                                            database='refseq_protein',
                                            sequence=sequence,
                                            entrez_query='txid117571 [ORGN]',  # boney fish & tetrapods
                                            hitlist_size=5000,
                                            )
        self.blast_record = self._blast_seq(sequence)
        self.results: pd.DataFrame = self._generate_df()

    def _blast_seq(self, sequence: str) -> Record.Blast:
        blast: io.StringIO = NCBIWWW.qblast(program='blastp',
                                            database='refseq_protein',
                                            sequence=sequence,
                                            entrez_query='txid117571 [ORGN]',  # boney fish & tetrapods
                                            hitlist_size=5000,
                                            )
        return list(NCBIXML.parse(blast))[0]

    def _generate_df(self) -> pd.DataFrame:
        results = pd.DataFrame(dict(title=map(operator.attrgetter('title'), self.blast_record.alignments),
                                    hsp=[al.hsps[0] for al in self.blast_record.alignments],
                                    ))
        results['ids'] = results.title.str.findall(r'\|(\w+_\d+\.\d+)\|(.*?)\[(.*?)\]')
        results['species'] = results.title.str.findall(r'\[(.*?)\]')
        results['species'] = results.title.str.findall(r'\[(.*?)\]')
        # assert sum(results['species'].apply(len) != results['ids'].apply(len)) == 0, 'There is a dogdy entry'
        parts = results.title.str.findall(r'\|(\w+_\d+\.\d+)\|\s?(.*?)\s?\[(.*?)\]')
        results['ids'] = parts.apply(functools.partial(map, operator.itemgetter(0))).apply(list)
        results['names'] = parts.apply(functools.partial(map, operator.itemgetter(1))).apply(list)
        results['specieses'] = parts.apply(functools.partial(map, operator.itemgetter(2))).apply(list)
        results['species'] = results.specieses.apply(operator.itemgetter(0))
        results['evalue'] = results.hsp.apply(operator.attrgetter('expect'))
        results['score'] = results.hsp.apply(operator.attrgetter('score'))
        results['matched_seq'] = results.hsp.apply(operator.attrgetter('sbjct')).str.replace('-', '')
        results['novel_species'] = results.specieses.apply(Novelty())
        results['fraction_identical'] = results.hsp.apply(operator.attrgetter('identities')) / results.hsp.apply(
            operator.attrgetter('align_length'))
        return results

    @classmethod
    def from_uniprot(cls, uniprot: str) -> Blaster:
        fasta = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.fasta').text
        seq = ''.join(fasta.split('\n')[1:])
        return cls(seq)

    clean = staticmethod(clean)
