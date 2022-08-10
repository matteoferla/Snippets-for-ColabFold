"""
Main function: align
"""

import pymol2, os, re
from matplotlib import cm
import numpy as np
from typing import Optional
from IPython.display import display, FileLink


def align(jobname, folder='.', bg_color: Optional[str] = None, **settings) -> None:
    """
    Given a jobname save (and display a link to) an pretty PyMol session.
    Chain A will be a shade of grey —ColabFold multimer does not do these
    Chain B onwards will be a shade of Greens, Blues, Reds, Purples, Oranges
    Model is model number. For binding affinity see ``pyrosetta-help`` package

    :param jobname: This includes the hash, for a path use folder.
    :param folder: The session is save in the current folder, whereas the data of jobname is from here.
    :param bg_color: optional color for the background. e.g. 'white'
    :param settings: any key value is passed to ``cmd.set()``
    :return:
    """
    modelnames = []
    with pymol2.PyMOL() as pymol:
        top_ranked = None
        model_names = []
        for filename in os.listdir(folder):
            if os.path.splitext(filename)[1] != '.pdb':
                continue
            match = re.match(jobname + r'_unrelaxed_rank_(?P<rank>\d)_model_(?P<model>\d+).pdb', filename)
            if not match:
                continue
            modelname = f'model{match.group("model")}'
            modelnames.append(modelname)
            if int(match.group('rank')) == 1:
                top_ranked = modelname
            pymol.cmd.load(os.path.join(folder, filename), modelname)
        num_models = len(modelnames)
        pymol.cmd.color('white', 'element C')
        print(f'job {jobname} has {num_models} models')
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
        for modelname in map('model{0}'.format, range(1, 6)):
            if modelname == top_ranked:
                continue
            pymol.cmd.align(f'{modelname} and chain B', f'{top_ranked} and chain B')
            for letter in subcolor_maps:
                pymol.cmd.show('sticks', f'{modelname} and not chain {letter} and ' +
                               f'byresi ({modelname} and chain {letter}) around 4')
                if letter == 'B':
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
