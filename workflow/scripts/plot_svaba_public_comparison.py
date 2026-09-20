#!/usr/bin/env python3
"""Plot measured truth-allele recall from the audited public comparison receipt."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from build_structural_exclusions import fingerprint, save, stamp


def plot(receipt, output, provenance):
    s=json.loads(receipt.read_text());overlap=s['truth_overlap_with_DELLY'];base=s['by_class']['baseline']
    if s['identical_truth_universe']!=802 or sum(x['total'] for x in overlap.values())!=802:
        raise ValueError('public truth accounting is incomplete')
    colors=['#57646d','#087f8c','#a9bfcc'];labels=['DELLY','SvABA','Either caller*']
    fig,axes=plt.subplots(1,2,figsize=(9.2,4.4),sharex=True,sharey=True)
    data=[]
    for ax,kind,title in zip(axes,['DEL','INS'],['Deletions','Insertions']):
        group='DEL' if kind=='DEL' else 'INS_like';denom=base[group+':TP']+base[group+':FN']
        both=overlap['matched_by_both']['classes'].get(kind,0)
        old=both+overlap['DELLY_only']['classes'].get(kind,0)
        new=both+overlap['SvABA_only']['classes'].get(kind,0)
        values=[old,new,old+new-both];rates=[100*x/denom for x in values]
        bars=ax.barh([2,1,0],rates,color=colors,height=.57)
        bars[-1].set_hatch('///');bars[-1].set_edgecolor('#718895')
        for y,n,rate in zip([2,1,0],values,rates):
            ax.text(rate+2,y,f'{n}/{denom}  ({rate:.1f}%)',va='center',fontsize=10)
        ax.set_title(title,fontsize=13,loc='left',pad=15,fontweight='bold')
        ax.set_yticks([2,1,0],labels);ax.set_xlim(0,100);ax.set_ylim(-.65,2.65)
        ax.set_xticks([0,25,50,75,100]);ax.set_xticklabels(['0%','25%','50%','75%','100%'])
        ax.set_xlabel('Known truth alleles matched',labelpad=10)
        ax.set_axisbelow(True);ax.grid(axis='x',color='#e4e8eb',linewidth=.8)
        ax.tick_params(axis='both',length=0,labelsize=10)
        for spine in ax.spines.values():spine.set_visible(False)
        data.append({'type':kind,'truth_alleles':denom,'matched':dict(zip(labels,values))})
    fig.suptitle('Complementary detection, substantial remaining gaps',x=.025,y=.985,ha='left',fontsize=16,fontweight='bold')
    fig.text(.025,.875,'Public HG002 chromosome 20 · GIAB v5.0q · fixed Truvari/MAFFT comparison',fontsize=10,color='#4c555c')
    fig.text(.025,.105,'* Matched by either in separate analyses; no merged-call precision has been measured.',fontsize=9)
    fig.text(.025,.055,'Single prealigned chromosome. This does not establish patient sensitivity or genotype accuracy.',fontsize=9,color='#4c555c')
    fig.subplots_adjust(left=.145,right=.965,top=.74,bottom=.265,wspace=.24)
    output.parent.mkdir(parents=True,exist_ok=True);artifacts=[]
    for extension in ('.png','.svg'):
        path=output.with_suffix(extension)
        if path.exists():raise FileExistsError(path)
        fig.savefig(path,dpi=180,facecolor='white');artifacts.append(fingerprint(path))
    plt.close(fig)
    result={'created_at':stamp(),'input':fingerprint(receipt),'code':fingerprint(Path(__file__)),
            'plotted_data':data,'artifacts':artifacts,'matplotlib_version':matplotlib.__version__}
    save(provenance,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('receipt','output','provenance'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();print(json.dumps(plot(a.receipt,a.output,a.provenance)))
