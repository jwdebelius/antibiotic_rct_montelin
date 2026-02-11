import numpy as np
import pandas as pd

colormap = {
    ('c__Bacteroidia', 'f__Bacteroidaceae'): '#084594',
    ('c__Bacteroidia', 'f__Rikenellaceae'): '#2171b5',
    ('c__Bacteroidia', 'f__Marinifilaceae'): '#4292c6',
    ('c__Bacteroidia', 'c__Bacteroidia'): '#9ecae1',
    ('c__Bacilli', 'f__Erysipelotrichaceae'): '#d94801',
    ('c__Bacilli', 'f__Erysipelatoclostridiaceae'): '#f16913',
    ('c__Bacilli', 'f__CAG-1000'): '#fd8d3c',
    ('c__Bacilli', 'c__Bacilli'): '#fdae6b',
    ('c__Clostridia', 'f__Lachnospiraceae'): '#99000d',
    ('c__Clostridia', 'f__Oscillospiraceae'): '#cb181d',
    ('c__Clostridia', 'f__Ruminococcaceae'): '#ef3b2c',
    ('c__Clostridia', 'f__Acutalibacteraceae'): '#fb6a4a',
    ('c__Clostridia', 'f__CAG-272'): '#fc9272',
    ('c__Clostridia', 'c__Clostridia'): '#fc9272',
    ('c__Negativicutes', 'f__Veillonellaceae'): '#6a51a3',
    ('c__Negativicutes', 'c__Negativicutes'): '#bcbddc',
    ('c__Other', 'f__Eggerthellaceae'): '#525252',
    ('c__Other', 'f__Desulfovibrionaceae'): '#737373',
    ('c__Other', 'c__Other'): '#bdbdbd',
    ('c__Gammaproteobacteria', 'f__Pasteurellaceae'): '#005a32',
    ('c__Gammaproteobacteria', 'f__Burkholderiaceae'): '#238b45',
    ('c__Gammaproteobacteria', 'f__Pseudomonadaceae'): '#41ab5d',
    ('c__Gammaproteobacteria', 'c__Gammaproteobacteria'): '#a1d99b',
}
color_order = [
    [
        ('c__Clostridia', 'f__Lachnospiraceae'),
         ('c__Clostridia', 'f__Oscillospiraceae'),
         ('c__Clostridia', 'f__Ruminococcaceae'),
         ('c__Clostridia', 'f__Acutalibacteraceae'),
         ('c__Clostridia', 'c__Clostridia'),
    ],
    [
        ('c__Negativicutes', 'f__Veillonellaceae'),
        ('c__Negativicutes', 'c__Negativicutes'),
    ],
    [
        ('c__Bacilli', 'f__Erysipelotrichaceae'),
        ('c__Bacilli', 'f__Erysipelatoclostridiaceae'),
        ('c__Bacilli', 'f__CAG-1000'),
        ('c__Bacilli', 'c__Bacilli'),
    ],
    [
        ('c__Bacteroidia', 'f__Bacteroidaceae'),
        ('c__Bacteroidia', 'f__Rikenellaceae'),
        ('c__Bacteroidia', 'f__Marinifilaceae'),
        ('c__Bacteroidia', 'c__Bacteroidia'),
    ],
    [
        ('c__Gammaproteobacteria', 'f__Pasteurellaceae'),
        ('c__Gammaproteobacteria', 'f__Burkholderiaceae'),
        ('c__Gammaproteobacteria', 'f__Pseudomonadaceae'),
        ('c__Gammaproteobacteria', 'c__Gammaproteobacteria'),
    ],
    [
        ('c__Other', 'f__Eggerthellaceae'),
        ('c__Other', 'f__Desulfovibrionaceae'),
        ('c__Other', 'c__Other'),
    ]
]
cmap_fam = pd.Series(colormap)
cmap_fam.reset_index(level=0, drop=True, inplace=True)

cmap2 = pd.Series(colormap).to_frame(name='color')
cmap2.index.set_names(['class', 'family'], inplace=True)
cmap2.reset_index(level='class', inplace=True)
cmap2 = cmap2.loc[np.hstack([[v[-1] for v in c] for c in color_order])]

cmap2['x'] = (~cmap2['class'].isin(['c__Clostridia', 'c__Negativicutes', 'c__Bacilli'])) * 1
cmap2['y'] = 1 + (~cmap2['class'].duplicated()) * 1.5
cmap2['y'] = cmap2.groupby('x')['y'].cumsum()
cmap2['y'] = cmap2['y'] - (cmap2['y'].min() - 1)

class_color = {
    'c__Clostridia': '#de2d26',
    'c__Bacteroidia': '#3182bd',
    'c__Gammaproteobacteria': '#31a354',
    'c__Negativicutes': '#756bb1',
    'c__Bacilli': '#e6550d',
    'c__Other': '#636363',
}
class_order = ['c__Clostridia', 'c__Negativicutes', 'c__Bacilli',
             'c__Bacteroidia', 'c__Gammaproteobacteria', 'c__Other']