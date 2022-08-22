"""
Main function: align
"""

import warnings
from typing import Optional

import numpy as np
import os
import pymol2
import re
from IPython.display import display, FileLink
from matplotlib import cm


def align(jobname: str,
          folder='.',
          align_chain='B',
          bg_color: Optional[str] = None, **settings):
    modelnames = []
    with pymol2.PyMOL() as pymol:
        top_ranked = None
        for filename in os.listdir(folder):
            if os.path.splitext(filename)[1] != '.pdb':
                continue
            match = re.match(jobname + r'_unrelaxed_rank_(?P<rank>\d)_model_(?P<model>\d+).pdb', filename)
            match_unfinished = re.match(jobname + r'_unrelaxed_model_(?P<model>\d+).pdb', filename)
            if match:
                modelname = f'model{match.group("model")}'
                if int(match.group('rank')) == 1:
                    top_ranked = modelname
            elif match_unfinished:
                warnings.warn('The analyses are not complete!')
                modelname = f'model{match_unfinished.group("model")}'
                top_ranked = modelname
            else:
                continue
            modelnames.append(modelname)
            pymol.cmd.load(os.path.join(folder, filename), modelname)
        # colour once found all
        num_models = len(modelnames)
        print(pymol.cmd.get_names())
        assert num_models > 0, f'there were no models for {jobname}'
        # technically the selection should return 5 atoms:
        assert pymol.cmd.select(f'chain {align_chain} and resi 1 and name CA'), f'There is not chain {align_chain}'
        pymol.cmd.color('white', 'element C')
        print(f'job {jobname} has {num_models} models')
        # 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds'
        # https://matplotlib.org/stable/gallery/color/colormap_reference.html
        color_maps = {
            'A': map(cm.get_cmap('Greys'), np.arange(0, 1, 1 / num_models)),  # colab multimer does not use A.
            'B': map(cm.get_cmap('Greens'), np.arange(0, 1, 1 / num_models)),
            'C': map(cm.get_cmap('Blues'), np.arange(0, 1, 1 / num_models)),
            'D': map(cm.get_cmap('Reds'), np.arange(0, 1, 1 / num_models)),
            'E': map(cm.get_cmap('Purples'), np.arange(0, 1, 1 / num_models)),
            'F': map(cm.get_cmap('Oranges'), np.arange(0, 1, 1 / num_models)),
        }
        n_chains = pymol.cmd.select('name CA and resi 1 and model1')
        subcolor_maps = {chr(66 + i): color_maps[chr(66 + i)] for i in range(n_chains)}
        for modelname in modelnames:
            for chain, color_map in subcolor_maps.items():
                pymol.cmd.color(cm.colors.to_hex(next(color_map)).replace('#', '0x'),
                                f'{modelname} and chain {chain} and element C')
        for modelname in map('model{0}'.format, range(1, num_models + 1)):
            if modelname == top_ranked or num_models == 1:  # or status_unfinished:
                continue
            pymol.cmd.align(f'{modelname} and chain {align_chain}', f'{top_ranked} and chain {align_chain}')
            for letter in subcolor_maps:
                pymol.cmd.show('sticks', f'{modelname} and not chain {letter} and ' +
                               f'byresi ({modelname} and chain {letter}) around 4')
                if letter == align_chain:
                    continue
                if pymol.cmd.select(f'{top_ranked} and chain {letter}') == 0 or pymol.cmd.select(
                        f'{modelname} and chain {letter}') == 0:
                    continue
                rmsd = pymol.cmd.rms_cur(f'{top_ranked} and chain {letter}', f'{modelname} and chain {letter}')
                print(f'Chain {letter} of {modelname} vs. {top_ranked}: {rmsd:.1f}')
        if bg_color:
            pymol.cmd.bg_color(bg_color)
        for setting, value in settings.items():
            pymol.cmd.set(setting, value)
        pymol.cmd.save(f'{jobname}.pse')

    # Green blue red
    display(FileLink(f'{jobname}.pse'))
