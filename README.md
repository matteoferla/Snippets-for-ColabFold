# Snippets-for-ColabFold
A collections of Pythonic snippets for working with ColabFold.

## Note
I am a biochemist not a geneticist.
I write in Python and do not use R. Please do not email me asking for a R snippet.

## Second note
These may be handy:

A Gist to convert a colab notebook to a Jupyter notebook:
[Gist: matteoferla/colab2jupyter.py](https://gist.github.com/matteoferla/15b483ab8f0c78293606dad91a360f9e)

A pip-module to import clean a given variable from a Gist, which features below.
[pypi: gist-import](https://pypi.org/project/gist-import/)
> pip install gist-import

## Files herein

### AnalyseA3M

Don't blindly blame your MSA, Analyse what is in a a3m. [AnalyseA3M](analyse_a3m.py).
This fetches Uniprot details of the entries and returns a pandas DataFrame for easy analysis.

```python
from gist_import import GistImporter
gi = GistImporter.from_github('https://github.com/matteoferla/Snippets-for-ColabFold/blob/main/analyse_a3m.py')
AnalyseA3M = gi['AnalyseA3M']

a3m = AnalyseA3M('VDAC2_BAK_tBID_2c974.a3m')
a3m.load_uniprot_data()
omnia = a3m.to_df()
boney = a3m.to_boney_subset()
print(f'{len(boney)} out of {len(omnia)} are tetrapods & boney-fish')
a3m.display_name_tallies(a3m.to_boney_subset())
boney.sample(5)
```
And do whatever filtering:

```python
# names are messy...
cleaned = boney.name_B.str.lower()\
                      .str.replace(r'( ?\d+)$','')\
                      .str.replace(r' proteinous','')\
                      .str.replace(r' protein','')\
                      .str.replace(r'bcl 2','bcl2')\
                      .str.strip()

# some homologues have question marks...
filtro = cleaned.isin(['bak', 'bcl2 antagonist/killer',
                           'bax regulator', 
                           #'bcl2', 'apoptosis regulator bcl 2',  # not sure.
                          # 'apoptosis regulator bcl',
                          # 'bcl domain containing'
                         ])
subsetted: pd.DataFrame = boney.loc[filtro & (boney.name_C == 'BH3 interacting domain death agonist')]
a3m.dump_a3m(subsetted, 'VDAC2_BAK_tBID_filtered.a3m')
```

### PyMOL alignment
Make pretty multimodel PyMOL alignment: [pymol_assemble](pymol_assemble.py)
```python
align(jobname='VDAC2_BAK_tBID_3693c', bg_color='white', use_shaders=0, ray_trace_mode=3)
```

### Chop up
Chop up an alphafold as a series of template with pLDDT > 70%: [chop](chop.py)

```python
folder = 'templates'

split_n_save('P45880', folder)
```
## PyRosetta

PyRosetta notebooks for working with AF2 models are in:

[GitHub: matteoferla/pyrosetta_help](https://github.com/matteoferla/pyrosetta_help)

These include scoring interface, phosphorylating, stretching the termini etc.

It also includes porting ligands over, but a way better version has been release (different group)
called [AlphaFill](https://alphafill.eu/)

## See also

Uniprot ported gnomAD variants:

[GitHub: matteoferla/Uniprot-ported-gnomADs](https://github.com/matteoferla/Uniprot-ported-gnomADs)

A series of PDBs for a blog post I will one day finish:

[GitHub: GitHub: matteoferla/autophagic-cell-death-complex-models](https://github.com/matteoferla/autophagic-cell-death-complex-models)


## Cluster Footnote

A collection of snippets for working with the University of Oxford's rescomp cluster (a SGE job scheduler system):
[GitHub: matteoferla/rescomp-tests](https://github.com/matteoferla/rescomp-tests)

But in terms of ColabFold on the cluster,
you need the internet for the MMSeq2 step (or make your own MSA). However, if you request zero models where there is 
internet it will generate a A3M. This can be used as a custom job on a node with no internet.
Note there is a hackish fix for now required — see below.

All my work in Rescomp is on a jupyter notebook forwarded from a private node.
Do not run notebooks in the log-in node. Sure, the R folk do it all the time, but do be a better citizen.
Get an interactive job instead.

```python
%%rescomp jobname=xxxx queue=short.qc cores=6
# needs to have had import rescomp on previous cell.

# hashed:
jobname = 'xxxx'
print(jobname)
print('no templates')
msa_mode = 'custom'
a3m_file = f"{jobname}.a3m"
result_dir="."
num_recycles = 3 # [1,3,6,12,24,48]
num_models = 5 # 5
dpi = 200

# ---------------------------------------------
import sys

from colabfold.download import download_alphafold_params, default_data_dir
from colabfold.utils import setup_logging
from colabfold.batch import get_queries, run, set_model_type
from colabfold.colabfold import plot_protein
from pathlib import Path
import matplotlib.pyplot as plt

setup_logging(Path(".").joinpath("log.txt"))
queries_path=f"{jobname}.csv"
queries_path=a3m_file
queries, is_complex = get_queries(queries_path)
run(
        queries=queries,
        result_dir='results',
        use_templates=False,
        custom_template_path=None,
        use_amber=False,
        msa_mode='custom',    
        model_type= "AlphaFold2-multimer-v2",
        num_models=num_models,
        num_recycles=num_recycles,
        model_order=[1, 2, 3, 4, 5],
        is_complex=is_complex,
        data_dir=Path("../ColabFoldData"),
        keep_existing_results=False,
        recompile_padding=1.0,
        rank_by="auto",
        pair_mode="unpaired+paired",
        stop_at_score=float(100),
        #prediction_callback=prediction_callback,
        dpi=dpi
)
```
## Shut up

If you are runnng 3.8 or above (Colab is stuck on 3.6), jax will yell at you things like
`FutureWarning: jax.tree_flatten is deprecated`. Supress it thusly:

```python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
```

## Hackish fix

Custom MSA in ColabFold version x.x.x blows up when paired alignments are passed.
Here is a fix:

```python
# hack the sourcecode!
from colabfold import batch

with open(batch.__file__, 'r') as fh:
    code = fh.read()
    
with open(batch.__file__.replace('.py', '.bk.py'), 'w') as fh:
    fh.write(code)
          
with open(batch.__file__, 'w') as fh:
    fh.write(code.replace('if unpaired_msa is None:', 'if unpaired_msa is None or unpaired_msa[sequence_index] == "":'))
```
