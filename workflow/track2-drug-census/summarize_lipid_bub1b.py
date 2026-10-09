#!/usr/bin/env python3
"""Reproduce a descriptive BUB1B extraction from deposited GSE220856 data.

This does not recreate count-based differential expression, calculate a P value,
or establish protein or functional rescue. Inputs remain unmodified.
"""
import argparse,csv,gzip,hashlib,json,math,re,statistics
from pathlib import Path

MATRIX_SHA='4dc1efaa03af0bb7c9d767122106f519b3c1dd864387a02cd63361f4020d1e63'
SOFT_SHA='5ba36f5da83416c9b029011e29c4385ad8bb49fd18c636caae5b5d534dc27cde'
GENE='ENSG00000156970'
GROUPS=[('control','C'),('PA','PA'),('PA+OA','OA'),('PA+BMS','BMS')]

def read_checked(path,expected):
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=expected:
        raise ValueError('Source snapshot differs: '+path.name)
    return gzip.decompress(data).decode()

def summarize(matrix,soft):
    rows=[r for r in csv.DictReader(matrix.splitlines(),delimiter='\t') if r['gene_id']==GENE]
    if len(rows)!=1 or rows[0]['gene_name']!='BUB1B':
        raise ValueError('Expected exact BUB1B feature, not the BUB1B-PAK6 readthrough')
    row=rows[0];samples={}
    for block in re.split(r'(?m)^\^SAMPLE = ',soft)[1:]:
        lines=block.splitlines();gsm=lines[0]
        title=[s.partition(' = ')[2] for s in lines if s.startswith('!Sample_title = ')]
        treatment=[s.partition('treatment: ')[2] for s in lines if s.startswith('!Sample_characteristics_ch1 = treatment: ')]
        if len(title)!=1 or len(treatment)!=1 or title[0] in samples:
            raise ValueError('Ambiguous sample mapping')
        samples[title[0]]={'GSM':gsm,'treatment':treatment[0]}
    if len(samples)!=12:raise ValueError('Expected12 deposited RNA samples')
    groups={};mapping=[]
    for label,prefix in GROUPS:
        columns=[];values=[]
        for i in range(1,4):
            title=f'HepaRG cells, {label}-{i}'
            if title not in samples or samples[title]['treatment']!=label:
                raise ValueError('Treatment/sample mapping differs: '+title)
            column=prefix+str(i);value=float(row[column])
            if not math.isfinite(value) or value<0:raise ValueError('Invalid expression value')
            columns.append(column);values.append(value)
            mapping.append({'column':column,'GSM':samples[title]['GSM'],'title':title,'condition':label})
        groups[label]={'source_columns':columns,'values':values,'n_libraries':3,
                       'mean':statistics.mean(values),'sample_SD':statistics.stdev(values)}
    for group in groups.values():
        group['ratio_to_PA_mean']=group['mean']/groups['PA']['mean']
        group['ratio_to_control_mean']=group['mean']/groups['control']['mean']
    return {'dataset':'GSE220856','gene_id':GENE,'gene_name':'BUB1B',
            'units':'FPKM in filename/method description; one GEO format sentence says RPKM',
            'groups':groups,'sample_mapping':mapping,
            'source_sha256':{'GSE220856_gene_fpkm.txt.gz':MATRIX_SHA,'GSE220856_family.soft.gz':SOFT_SHA},
            'scope':'Descriptive extraction of deposited normalized values; no raw-count differential-expression reproduction, new P value or adjusted P value.',
            'limitations':['Three libraries per group, not independent donors.',
                           'PA and PA+OA have4h exposure; PA+BMS has24h BMS with concomitant4h PA.',
                           'No BMS-alone or EPA RNA samples are deposited in this series.',
                           'All three replicates retained, including BMS3.',
                           'No BUBR1 protein, chromosome fidelity, MVA or muscle rescue established.']}

def plot(result,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.4,5.2),dpi=170)
    colors=['#6b7280','#c5653b','#53947c','#336eab']
    labels=['Control','Palmitate','Palmitate\n+ oleic acid','Palmitate\n+ BMS-309403']
    for i,((group,_),color) in enumerate(zip(GROUPS,colors)):
        g=result['groups'][group]
        ax.scatter([i-.10,i,i+.10],g['values'],s=52,color=color,zorder=3)
        ax.hlines(g['mean'],i-.25,i+.25,color=color,lw=2.6)
    ax.set_xticks(range(4),labels);ax.set_ylabel('BUB1B expression (reported FPKM)')
    ax.set_ylim(0,23);ax.set_xlim(-.6,3.6)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True)
    ax.set_title('BUB1B RNA in a liver-cell injury model',loc='left',fontsize=15,pad=26)
    ax.text(0,1.03,'Deposited GSE220856 values • HepaRG cells',transform=ax.transAxes,color='#555',fontsize=10)
    fig.text(.12,.06,'Each point is one library; lines show means. No significance claim.\nBMS has a different exposure schedule; this is not protein or MVA rescue.',fontsize=9,color='#555')
    fig.subplots_adjust(left=.12,right=.98,bottom=.23,top=.82)
    fig.savefig(path,facecolor='white');plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True);p.add_argument('--soft',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--figure',type=Path)
    a=p.parse_args();result=summarize(read_checked(a.matrix,MATRIX_SHA),read_checked(a.soft,SOFT_SHA))
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    if a.figure:plot(result,a.figure)
    print(json.dumps({'gene':GENE,'libraries':12,'PA_BMS_over_PA_mean':result['groups']['PA+BMS']['ratio_to_PA_mean'],'inferential_statistics_computed':False},indent=2))

if __name__=='__main__':main()
