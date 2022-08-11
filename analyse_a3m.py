import itertools
import json
import operator
import re
from typing import Optional

import pandas as pd
import requests
from Bio import SeqIO
from Bio import pairwise2
from IPython.display import display, HTML


def retrieve_uniprot_data(uniprot: str) -> dict:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.json')
    response.raise_for_status()
    return response.json()


class AnalyseA3M:
    """
    Given a ColabFold alignment via MMseqs2,
    look up what the constituents are in order to debug with a complex failed to be predicted.

    This is primarily intended for human protein complexes, which are a lot less tidy than prokaryotic ones.

    ```python
    a3m = AnalyseA3M('VDAC2_BAK_tBID_2c974.a3m')
    # no uniprot data fetched yet
    omnia: pd.DataFrame = a3m.to_df()  # full dataframe of pairs
    # uniprot data has been fetched
    a3m.dump_uniprot_data() # calling `a3m.load_uniprot_data()` circumvebnts the retrieval
    boney: pd.DataFrame = a3m.to_boney_subset()
    print(f'{len(boney)} out of {len(omnia)} are tetrapods & boney-fish')
    # do some check etc. with the pd.DataFrame
    boney.sample(5) # etc.
    # shortcut:
    a3m.display_name_tallies(a3m.to_boney_subset())
    # save
    a3m.dump_a3m(boney.loc[boney.name_C == 'BH3 interacting domain death agonist'], 'VDAC2_BAK_tBID_filtered.a3m')
    ```

    """

    def __init__(self, alignmnet_filename: str):
        self.pairs = []
        self.seqs = []
        self.uniprot_data = {}
        self.df = pd.DataFrame()
        with open(alignmnet_filename, 'r') as fh:
            self._first_line = next(fh)
            for record in SeqIO.parse(fh, "fasta"):
                self.pairs.append(record.description.split())
                self.seqs.append(str(record.seq).replace('\n', ''))
        assert len(self.pairs), 'There are no sequences...'
        self.n_protein = len(self.pairs[0])
        self.lengths = list(map(int, self._first_line[1:].split('\t')[0].split(',')))
        # last paired
        self.last_pair = 0
        for pair in self.pairs[1:]:
            if len(pair) != self.n_protein:
                break
            self.last_pair += 1

    def _retrieve_uniprots(self):
        for uniprot in itertools.chain(*self.pairs[1:self.last_pair]):
            try:
                if uniprot[:3] == 'UPI':  # environmental
                    continue
                if uniprot in self.uniprot_data:
                    continue
                self.uniprot_data[uniprot] = retrieve_uniprot_data(uniprot)
            except requests.HTTPError as error:
                self.uniprot_data[uniprot] = None
                print(error)
        return self.uniprot_data

    def __getitem__(self, uniprot):
        return self.uniprot_data.get(uniprot, None)

    def load_uniprot_data(self, filename='uniprot_data.json'):
        """
        Loads a previously stored uniprot data retrieval.
        See ``.dump_uniprot_data``
        """
        with open(filename, 'r') as fh:
            self.uniprot_data.update(**json.load(fh))

    def dump_uniprot_data(self, filename='uniprot_data.json'):
        with open(filename, 'w') as fh:
            json.dump(self.uniprot_data, fh)

    def nested_item_getter(self, *keys):

        def getter(uniprot):
            sub = self[uniprot]
            if sub is None:
                return None
            for k in keys:
                if (isinstance(k, int) and len(sub) >= k) or (k in sub):
                    sub = sub[k]
                else:
                    return None
            return sub

        return getter

    def _calc_identity(self, seqs, i: int = -1):
        """
        A wrapper for ``Bio.pairwise2.identity_match``.
        which is a fun that returns 1 for an identity and 0 for a mismatch.
        (So the same as ``int(a == b)``). One day I might use a more sophisticated
        comparision metric.
        """
        if i == -1:
            begin = None
            end = None
        else:
            begin = sum(self.lengths[:i])
            end = begin + self.lengths[i] + 1
        r = slice(begin, end)
        # insertions in a3m are lowercase
        clean_seq = seqs.apply(lambda s: re.sub('[^A-Z-]', '', s)[r])
        # I want iterations per aa of interations per sequence
        matcher = pairwise2.identity_match()
        # dumb, but this get confusing quickly below...
        packed_matcher = lambda c: matcher(*c)  # noqa: E731 Guido doesn't like this, but I do
        seq_matcher = lambda pair: sum(map(packed_matcher, pair))  # noqa: E731 sorry Guido
        s_iter = clean_seq.apply(lambda s: zip(s, self.seqs[0][r]))
        return s_iter.apply(seq_matcher) / len(self.seqs[0][r].replace('-', ''))

    def to_df(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame and stores it in the attribute ``df``.
        Do note, the DataFrame is calculated only if ``df`` is empty, e.g.

        ```python
        a3m = AnalyseA3M('foo.a3m')
        _ = a3m.to_df()  # calculated for first time
        _ = a3m.to_df()  # not
        a3m.df = pd.DataFrame()
        _ = a3m.to_df()  # recalculated
        ```
        """
        if len(self.df) > 0:
            return self.df

        if not self.uniprot_data:
            self._retrieve_uniprots()

        to_tuple_series = lambda A, B: pd.Series(zip(A.values, B.values))   # noqa: E731
        combine = lambda A, B: to_tuple_series(A, B).apply(lambda ab: ab[0] if ab[0] else ab[1])  # noqa: E731
        # make table
        df = pd.DataFrame(dict(pair=self.pairs[1:self.last_pair + 1], seq=self.seqs[1:self.last_pair + 1]))
        df['identity'] = self._calc_identity(df.seq)
        for i in range(self.n_protein):
            letter = chr(i + 65)
            df[letter] = df.pair.apply(operator.itemgetter(i))
            df[f'identity_{letter}'] = self._calc_identity(df.seq, i)
            submitted = df[letter].apply(
                self.nested_item_getter('proteinDescription', 'submissionNames', 0, 'fullName', 'value'))
            recommended = df[letter].apply(
                self.nested_item_getter('proteinDescription', 'recommendedName', 'fullName', 'value'))
            df[f'name_{letter}'] = combine(recommended, submitted).fillna('Unnamed')
            # , regex=True does not work on .str.extract
            df[f'number_{letter}'] = df[f'name_{letter}'].fillna('').astype(str).str.extract(r'(\d+)')[0].astype(float)
        df['species'] = df.A.apply(self.nested_item_getter('organism', 'scientificName'))
        for i, taxon in enumerate(('domain', 'kingdom', 'superphylum', 'phylum', 'subphylum', 'order')):
            df[taxon] = df.A.apply(self.nested_item_getter('organism', 'lineage', i))
        self.df = df
        return self.df

    def clean_name(self, name_series: pd.Series, genename: Optional[str] = None) -> pd.Series:
        """
        Returns a lowercase series of names without words like isoform X3, homolog or protein.
        If a gene name is passed, then any variant with spaces or similar is removed. e.g. B.L.T. -> blt
        """
        if genename:
            pattern = re.compile(''.join(map(lambda c: f'{c}\W?', str(genename))))  # noqa: it is valid.
        else:
            pattern = re.compile('This is just fluff that will never match')
        return name_series.fillna('Unnamed') \
            .astype(str) \
            .str.lower() \
            .str.replace(r'isoform .*', '', regex=True) \
            .str.replace('like', '') \
            .str.replace('homologous', '') \
            .str.replace('homologue', '') \
            .str.replace('homolog', '') \
            .str.replace('-', ' ') \
            .str.replace(r'( ?\d+)$', '', regex=True) \
            .str.replace(r' proteinous', '') \
            .str.replace(r' protein', '') \
            .str.replace(pattern, '', regex=True) \
            .str.replace('  ', ' ') \
            .str.strip()

    def _check_df(self, df: Optional[pd.DataFrame] = None, wanted=()) -> pd.DataFrame:
        if df is not None:
            pass
        elif len(self.df):
            df = self.df
        else:
            df = self.to_df()
        assert len(df), 'Can only work w/ non-zero length dataframes'
        for colname in wanted:
            assert colname in df.columns, f'This is a nonstandard dataframe: no {colname}'
        return df

    def to_boney_subset(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        I am unable to remember the correct spelling of Euteleostomi
        also who calls them that?
        Boney animals are fancy-fish (not sharks, rays and lampreys) + tetrapods
        """
        df = self._check_df(df, wanted=('order',))
        euteleostomi = df.loc[df.order == 'Euteleostomi']
        return euteleostomi

    def to_nonboney_subset(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        The opposite of ``to_boney_subset``.
        Boney animals are fancy-fish (not sharks, rays and lampreys) + tetrapods
        """
        df = self._check_df(df, wanted=('order',))
        euteleostomi = df.loc[df.order != 'Euteleostomi']
        return euteleostomi

    def dump_a3m(self, df: pd.DataFrame, filename: str = 'filtered.a3m'):
        df = self._check_df(df, wanted=('seq', 'name_A'))
        """Save the a3m of the filtered pandas dataframe.

        Why is df not optional? Unlike other commands the ``df`` argument is not optional
        as it only makes sense if you specify a subset of ``self.df``."""
        with open(filename, 'w') as fh:
            fh.write(self._first_line.strip() + '\n')
            first_header = "\t".join(self.pairs[0])
            fh.write(f'>{first_header}\n{self.seqs[0]}\n')
            uniprot_cols = [str(chr(65 + i)) for i in range(self.n_protein)]
            for i, row in df.iterrows():
                header: str = '\t'.join([row[uniprot_col] for uniprot_col in uniprot_cols])
                fh.write(f'>{header}\n{row.seq}\n')

    def display_name_tallies(self, df: Optional[pd.DataFrame] = None):
        """Display (Jupyter notebook) the number of times
        a given protein name is found."""
        name_cols = [f'name_{chr(65 + i)}' for i in range(self.n_protein)]
        df: pd.DataFrame = self._check_df(df, wanted=name_cols)
        for i in range(self.n_protein):
            letter = chr(i + 65)
            display(HTML(f'<h3>Tally of protein #{i} ({letter})</h3>'))
            display(df[f'name_{letter}'].value_counts())
