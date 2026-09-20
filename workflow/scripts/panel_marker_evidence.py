"""Qualify constrained marker likelihoods using their original pileup evidence.

Thresholds are declared research choices, not validated clinical detection limits.
The original caller result is preserved; exclusion never rewrites it as reference.
"""


def assess_marker(expected_alleles,pileup,call,minimum_depth=10,maximum_depth=100,
                  nonpanel_minimum_reads=3,nonpanel_minimum_fraction=0.10):
    if minimum_depth<1 or maximum_depth<minimum_depth or nonpanel_minimum_reads<1 or not 0<nonpanel_minimum_fraction<=1:
        raise ValueError('invalid marker evidence thresholds')
    reasons=[];metrics={}
    if tuple(call['alleles'])!=tuple(expected_alleles):
        raise ValueError('constrained alleles differ from independently selected marker')
    sample=call['sample'];dp=sample.get('DP');pl=sample.get('PL')
    if dp is None or dp==0:reasons.append('no_usable_depth')
    elif type(dp) is not int or dp<0:raise ValueError('invalid genotype depth')
    elif not minimum_depth<=dp<=maximum_depth:reasons.append('depth_outside_declared_range')
    if (pl is None or len(pl)!=3 or any(type(x) is not int or x<0 for x in pl)):
        reasons.append('missing_or_invalid_diploid_likelihoods')
    elif min(pl)!=0 or len(set(pl))==1:
        reasons.append('non_normalized_or_uninformative_likelihoods')
    if pileup is None:
        reasons.append('no_original_pileup_record')
    else:
        pdp=pileup['sample'].get('DP');ad=pileup['sample'].get('AD')
        if dp!=pdp:raise ValueError('constrained and original pileup depth disagree')
        if ad is None or len(ad)!=len(pileup['alleles']) or any(type(x) is not int or x<0 for x in ad):
            reasons.append('missing_or_invalid_original_allele_depths')
        else:
            offpanel=0;unknown=0
            for allele,depth in zip(pileup['alleles'],ad):
                if allele in expected_alleles:continue
                if len(allele)==1 and allele in 'ACGT':offpanel+=depth
                elif depth:unknown+=depth
            metrics['nonpanel_nucleotide_depth']=offpanel
            metrics['unresolved_allele_depth']=unknown
            metrics['nonpanel_nucleotide_fraction']=offpanel/dp if dp else None
            if unknown:reasons.append('unresolved_original_allele_counts')
            if dp and offpanel>=nonpanel_minimum_reads and offpanel/dp>=nonpanel_minimum_fraction:
                reasons.append('substantial_third_allele_evidence')
    return {'usable_for_declared_ROH_analysis':not reasons,'exclusion_reasons':reasons,
            'original_pileup_metrics':metrics,'thresholds':{'minimum_depth':minimum_depth,
            'maximum_depth':maximum_depth,'nonpanel_minimum_reads':nonpanel_minimum_reads,
            'nonpanel_minimum_fraction':nonpanel_minimum_fraction},
            'scope':'Likelihood/read-count qualification only; excluded markers remain unassessed, not reference or a clinical negative.'}
