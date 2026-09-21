#!/usr/bin/env python3
"""Render the completed results; this script does not run simulations."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / 'population-results.json').read_text())
alleles = json.loads((ROOT / 'allele-inputs.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
                     'axes.spines.top':False, 'axes.spines.right':False,
                     'axes.titleweight':'bold', 'axes.titlesize':12,
                     'figure.facecolor':'white', 'savefig.facecolor':'white'})
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.subplots_adjust(left=.19, right=.97, top=.84, bottom=.22, hspace=.8, wspace=.85)
blue, red, gold = '#226b95', '#af4936', '#b47c26'

ax = axes[0,0]
keys = ['WT', 'N1002Lys', 'L737Ter', 'L737Ter_N1002Lys_cis']
labels = ['Wild type', 'N1002K alone', 'L737* alone', 'Both variants in cis']
for i, key in enumerate(keys):
    n = alleles['allele_products'][key]['protein_length_aa']
    ax.barh(i, n, color=blue if n==1050 else gold, height=.55)
    ax.text(n+18, i, str(n), va='center', fontsize=10)
ax.plot([1002,1002],[.69,1.31], color=red, lw=3)
ax.set(yticks=range(4), yticklabels=labels, xlim=(0,1200), xticks=[0,500,1000],
       xlabel='Predicted translated length (amino acids)')
ax.invert_yaxis()
ax.set_title('A  Sequence constrains rescue', loc='left', pad=15)
ax.text(0, -.38, 'Conditional on translation; RNA and protein\nsurvival remain unknown.',
        transform=ax.transAxes, fontsize=9, color='#555555', wrap=True)

ax=axes[0,1]
case=data['illustrative_trajectories']['renewing']['arms']
base, kill=case['untreated'],case['selective_aneuploid_killing']
values=[100*(kill['euploid_fraction']-base['euploid_fraction']),
        100*(kill['raw_differentiated_cell_time']/base['raw_differentiated_cell_time']-1),
        100*(kill['weighted_functional_cell_time']/base['weighted_functional_cell_time']-1)]
ax.barh(range(3),values,color=[blue,red,red],height=.55)
ax.axvline(0,color='#777777',lw=.8)
for i,v in enumerate(values):
    ax.text(v+(.6 if v>0 else -.6),i,f'{v:+.1f}'+(' pp' if i==0 else '%'),
            va='center',ha='left' if v>0 else 'right',fontsize=10)
ax.set(yticks=range(3),yticklabels=['Euploid fraction\n(final)','Differentiated\ncell-time','Weighted\nfunctional cell-time'],
       xlim=(-18,24),xticks=[-10,0,10,20],xlabel='Change versus untreated\n(units shown on bars)')
ax.invert_yaxis()
ax.set_title('B  Enrichment can hide cell loss',loc='left',pad=15)
ax.text(0,-.38,'Assumed selective aneuploid killing;\nrenewing example, 20 units.',
        transform=ax.transAxes,fontsize=9,color='#555555')

ax=axes[1,0]
conditions=[(0,0),(0,.25),(1,0),(1,.25)]
for horizon,color,offset in [(20,blue,-.13),(80,gold,.13)]:
    vals=[]
    for harm,function in conditions:
        r=next(r for r in data['rows'] if r['duration']==horizon
               and r['intervention']=='selective_senescent_removal'
               and all(r['scenario'][k]==v for k,v in dict(renewal=True,
                   aneuploid_division_ratio=1.8,error_probability=.04,
                   initial_senescence_fraction=.25,senescence_harm=harm,
                   senescent_function=function).items()))
        vals.append(100*(r['weighted_output_ratio']-1))
    ax.scatter(vals,[i+offset for i in range(4)],color=color,s=44,label=f'{horizon} units',zorder=3)
ax.axvline(0,color='#777777',lw=.8)
ax.set(yticks=range(4),yticklabels=['No harm / no function','No harm / some function',
        'Neighbor harm / no function','Neighbor harm / some function'],
       xlim=(-10,30),xlabel='Change in weighted functional cell-time (%)')
ax.invert_yaxis()
ax.legend(frameon=False,loc='upper right',fontsize=9)
ax.set_title('C  Senolysis can change sign',loc='left',pad=15)
ax.text(0,-.4,'Same two pulses; no senescent niche occupancy.\n"Some function" = 0.25; neighbor-harm coefficient = 1.',
        transform=ax.transAxes,fontsize=9,color='#555555')

ax=axes[1,1]
maximum=[]
for horizon in data['horizons']:
    rs=[r for r in data['rows'] if r['duration']==horizon and r['intervention']=='growth_slowing_only']
    maximum.append(100*(max(r['required_differentiated_function_multiplier'] for r in rs)-1))
ax.bar([0,1],maximum,color=[blue,gold],width=.55)
for i,v in enumerate(maximum):
    ax.text(i,v+.5,f'{v:.1f}%',ha='center',fontsize=11)
ax.set(xticks=[0,1],xticklabels=['20 units','80 units'],ylim=(0,16),
       ylabel='Extra function per differentiated cell (%)',xlabel='Observation horizon\n(arbitrary model units)')
ax.set_title('D  Slower growth has a cost',loc='left',pad=15)
ax.text(0,-.4,'Largest required gain after 35% slower division.\nNot an estimate of mTOR-inhibitor benefit or harm.',
        transform=ax.transAxes,fontsize=9,color='#555555')

fig.suptitle('MVA: what the genome and population models can resolve',x=.08,y=.97,ha='left',fontsize=18,weight='bold')
fig.text(.08,.92,'Sequence reconstruction plus assumption sensitivity — not predicted drug responses',fontsize=12,color='#555555')
fig.text(.08,.045,'Population rates, initial cell composition and intervention effects are assumed. No patient days or clinical effect sizes.\nPanels B–C use a declared illustrative setting; panel D shows extrema in the tested grid. Full definitions and limitations: README.md.',
         fontsize=10,color='#555555')
fig.savefig(ROOT/'model-findings.png',dpi=180)
fig.savefig(ROOT/'model-findings.pdf')
print('Rendered model-findings.png and model-findings.pdf from retained results.')
