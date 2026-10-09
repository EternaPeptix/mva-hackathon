#!/usr/bin/env python3
"""Replay all 80 deposited records; this cannot recover the published n=73 set.

Requires NumPy and SciPy. Official CSV:
https://ndownloader.figshare.com/files/39372206
Dataset DOI: 10.6084/m9.figshare.22081895.v4

The seven efficacy exclusions are unidentified, and the deposit differs from
published subgroup summaries. All results are nominal exploratory calculations,
not a reproduced trial primary analysis, new experiment or clinical advice.
"""
from pathlib import Path
import argparse, csv, hashlib, json, sys
import numpy as np
import scipy
from scipy import stats

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('csv', type=Path, help='Unmodified official raw-data.csv')
parser.add_argument('output', type=Path, help='New output JSON path; existing files are not replaced')
args = parser.parse_args()
source = args.csv.resolve()
expected = '66d2e92ec3a8d7834dc9a869223117b12bdf50549619f30fe71df320392cbcd0'
actual = hashlib.sha256(source.read_bytes()).hexdigest()
if actual != expected:
    parser.error('Source SHA-256 differs from the reviewed public CSV; no analysis performed.')
if args.output.exists():
    parser.error('Output already exists; choose a new path.')
with source.open(encoding='utf-8-sig', newline='') as handle:
    rows = list(csv.DictReader(handle))
if len(rows) != 80 or any(sum(r['Group'] == g and r['Sex'] == s for r in rows) != 20
                        for g in ['KSM-66', 'Placebo'] for s in ['Male', 'Female']):
    parser.error('Expected all 80 records, with 20 per sex per arm.')
out = {'status':'INDEPENDENT_RECOMPUTATION_NOT_PUBLISHED_N73_ANALYSIS',
       'source':str(source), 'source_sha256':actual,
       'runtime':sys.executable, 'scipy_version':scipy.__version__, 'numpy_version':np.__version__,
       'analysis_set':'All 80 publicly deposited complete records; no exclusions or imputation. The seven compliance exclusions are not identified by the public columns.',
       'tests':'Two-sided independent pooled-variance Student t test and Welch sensitivity test, 95% nominal CIs; no multiplicity correction.',
       'percentage_definition':'100*(Day56-Day0)/Day0, recomputed from unrounded deposited measurement values; percentage-point contrast.',
       'outcomes':[], 'baseline_means':[]}

def contrast(x,y):
    n,m=len(x),len(y)
    v,w=np.var(x,ddof=1),np.var(y,ddof=1)
    delta=float(np.mean(x)-np.mean(y))
    sp=np.sqrt(((n-1)*v+(m-1)*w)/(n+m-2))
    pooled=stats.ttest_ind(x,y,equal_var=True)
    welch=stats.ttest_ind(x,y,equal_var=False)
    def result(test,se):
        width=stats.t.ppf(.975,test.df)*se
        return {'p':float(test.pvalue),'t':float(test.statistic),'df':float(test.df),'ci95':[float(delta-width),float(delta+width)]}
    return {'n_ag':n,'n_placebo':m,'ag_mean':float(np.mean(x)),'ag_sd':float(np.std(x,ddof=1)),
            'placebo_mean':float(np.mean(y)),'placebo_sd':float(np.std(y,ddof=1)),
            'difference_ag_minus_placebo':delta,'cohen_d_pooled':float(delta/sp),
            'student':result(pooled,sp*np.sqrt(1/n+1/m)),
            'welch':result(welch,np.sqrt(v/n+w/m))}

for name in ['1RM Chest Press','1RM Leg Press','Max. VO2','Mid-arm circumference','Thigh circumference','Chest circumference','Body fat %']:
    for sex in ['All','Male','Female']:
        group_rows=[[r for r in rows if r['Group']==group and (sex=='All' or r['Sex']==sex)] for group in ['KSM-66','Placebo']]
        before=[np.array([float(r[name+' (Day 0)']) for r in rr]) for rr in group_rows]
        after=[np.array([float(r[name+' (Day 56)']) for r in rr]) for rr in group_rows]
        out['baseline_means'].append({'outcome':name,'sex':sex,'ag':float(np.mean(before[0])),'placebo':float(np.mean(before[1]))})
        for scale in ['absolute_change','percentage_change']:
            change=[a-b if scale=='absolute_change' else 100*(a-b)/b for a,b in zip(after,before)]
            out['outcomes'].append({'outcome':name,'sex':sex,'scale':scale,**contrast(*change)})

# Reproduce the stated covariate structure as a labeled exploratory fit to all80,
# not a claim to reproduce Table4's unpublished n73 subset or software settings.
X=np.array([[1,1 if r['Group']=='KSM-66' else 0,float(r['Gender']),float(r['BMI (Day 0)']),float(r['Chest circumference (Day 0)'])] for r in rows])
out['all80_exploratory_ancova']=[]
for name in ['1RM Chest Press','1RM Leg Press','Max. VO2']:
    before=np.array([float(r[name+' (Day 0)']) for r in rows]); after=np.array([float(r[name+' (Day 56)']) for r in rows])
    for scale in ['absolute_change','percentage_change']:
        y=after-before if scale=='absolute_change' else 100*(after-before)/before
        beta=np.linalg.lstsq(X,y,rcond=None)[0]; df=len(y)-np.linalg.matrix_rank(X)
        cov=np.linalg.inv(X.T@X)*(sum((y-X@beta)**2)/df);se=np.sqrt(cov[1,1]);t=beta[1]/se
        width=stats.t.ppf(.975,df)*se
        out['all80_exploratory_ancova'].append({'outcome':name,'scale':scale,'formula':'change ~ treatment + numeric sex + baseline BMI + baseline chest circumference',
          'n':len(y),'df':int(df),'treatment_difference':float(beta[1]),'ci95':[float(beta[1]-width),float(beta[1]+width)],'p':float(2*stats.t.sf(abs(t),df))})
with args.output.open('x') as handle:
    json.dump(out, handle, indent=2, allow_nan=False)
    handle.write('\n')
print(f"Wrote {len(out['outcomes'])} exploratory contrasts for all 80 records; not the published n=73 analysis.")
