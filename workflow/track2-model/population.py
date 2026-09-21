#!/usr/bin/env python3
"""Genome-conditioned sensitivity analysis, not a fitted model of this child.

Rates are per arbitrary model-time unit. Genotype and karyotype are distinct.
No parameter set represents measured efficacy of a named drug.
"""
import argparse
import hashlib
import itertools
import json
import math
import platform
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class Scenario:
    renewal: bool
    aneuploid_division_ratio: float
    error_probability: float
    initial_senescence_fraction: float
    senescence_harm: float
    senescent_function: float
    senescent_niche_weight: float = 0.0
    aneuploid_function: float = 0.5
    senolytic_offtarget_fraction: float = 0.02
    capacity: float = 2000.0


INTERVENTIONS = (
    'untreated', 'functional_error_reduction', 'stress_protection',
    'growth_slowing_only', 'growth_slowing_and_less_senescence',
    'selective_senescent_removal', 'senescent_removal_with_offtarget_loss',
    'selective_aneuploid_killing',
)
HORIZONS = (20.0, 80.0)


def competition_regimes():
    # Untreated progenitor losses include differentiation, senescence and death.
    # At a hypothetical E-maintained equilibrium, A grows if r exceeds this
    # threshold. Mutation prevents a strictly E-only equilibrium when mu > 0.
    rows = []
    for mu in (0.005, 0.02, 0.04):
        threshold = (0.10 / 0.065) * (1.0 - 2.0 * mu)
        assert 0.7 < threshold < 1.8
        rows.append({'error_probability':mu, 'division_ratio_boundary':threshold,
                     'lower_ratio':0.7, 'upper_ratio':1.8})
    return {'scope':'Baseline progenitor selection boundary, conditional on the stated loss rates and a positive crowded equilibrium. Not measured biological fitness.',
            'formula':'r_boundary = (0.10 / 0.065) * (1 - 2*mu)', 'rows':rows}


def initial_state(s):
    # 4,000 starting cells; 20% abnormal in each compartment is an assumption.
    n = 4000.0
    sen = n * s.initial_senescence_fraction
    proliferating = min(800.0, n - sen) if s.renewal else 0.0
    differentiated = n - sen - proliferating
    return [0.8 * proliferating, 0.2 * proliferating,
            0.8 * differentiated, 0.2 * differentiated,
            0.8 * sen, 0.2 * sen, 0.0, 0.0, 0.0, 0.0]


def functional_components(y, s):
    pe, pa, de, da, se, sa = y[:6]
    # Raw counts are reported alongside this assumed function model.
    denominator = 1.0 + s.senescence_harm * (se + sa) / s.capacity
    return ((de + s.aneuploid_function * da) / denominator,
            s.senescent_function * (se + sa) / denominator)


def simulate(s, intervention='untreated', *, dt=0.05, duration=20.0,
             error_factor=0.5, retain_trajectory=False):
    if intervention not in INTERVENTIONS:
        raise ValueError(intervention)
    if not 0 <= error_factor <= 1 or not 0 <= s.error_probability <= 1:
        raise ValueError('invalid probability')
    steps = round(duration / dt)
    if not math.isclose(steps * dt, duration, abs_tol=1e-10):
        raise ValueError('duration must be an integer number of steps')
    y = initial_state(s)
    n0 = sum(y[:6])
    integral_d = integral_s = integral_raw_d = 0.0
    max_accounting_error = 0.0
    max_abnormal = y[1] + y[3] + y[5]
    history = []
    for step in range(steps):
        t = step * dt
        pe, pa, de, da, se, sa, dead, divisions, errors, transfers = y
        b = 0.30 if s.renewal else 0.0
        differentiation = 0.05 if s.renewal else 0.0
        p_sen_e, p_sen_a = 0.01, 0.03
        d_sen_e, d_sen_a = 0.002, 0.005
        p_die_e, p_die_a = 0.005, 0.02
        d_die_e, d_die_a = 0.005, 0.01
        s_die = 0.001
        mu = s.error_probability
        off_target = 0.0
        if intervention == 'functional_error_reduction':
            mu *= error_factor
        elif intervention == 'stress_protection':
            # Hypothetical protection of already abnormal cells, not error repair.
            p_sen_a *= 0.5
            d_sen_a *= 0.5
            p_die_a *= 0.5
            d_die_a *= 0.5
        elif intervention in ('growth_slowing_only', 'growth_slowing_and_less_senescence'):
            b *= 0.65
            if intervention == 'growth_slowing_and_less_senescence':
                p_sen_e *= 0.5
                p_sen_a *= 0.5
                d_sen_e *= 0.5
                d_sen_a *= 0.5
        elif intervention in ('selective_senescent_removal', 'senescent_removal_with_offtarget_loss'):
            # Two assumed pulses, not a clinical schedule.
            pulse = (5.0 <= t < 6.0) or (15.0 <= t < 16.0)
            if pulse:
                s_die += math.log(2)
                if intervention == 'senescent_removal_with_offtarget_loss':
                    off_target = -math.log1p(-s.senolytic_offtarget_fraction)
        elif intervention == 'selective_aneuploid_killing':
            p_die_a += 0.1
            d_die_a += 0.1
        p_die_e += off_target
        p_die_a += off_target
        d_die_e += off_target
        d_die_a += off_target
        occupied = pe + pa + 0.1 * (de + da) + s.senescent_niche_weight * (se + sa)
        gate = max(0.0, 1.0 - occupied / s.capacity)
        be = b * gate * pe
        ba = b * s.aneuploid_division_ratio * gate * pa
        se_in = p_sen_e * pe + d_sen_e * de
        sa_in = p_sen_a * pa + d_sen_a * da
        deaths = p_die_e * pe + p_die_a * pa + d_die_e * de + d_die_a * da + s_die * (se + sa)
        dy = [be * (1.0 - 2.0 * mu) - (differentiation + p_sen_e + p_die_e) * pe,
              2.0 * mu * be + ba - (differentiation + p_sen_a + p_die_a) * pa,
              differentiation * pe - (d_sen_e + d_die_e) * de,
              differentiation * pa - (d_sen_a + d_die_a) * da,
              se_in - s_die * se, sa_in - s_die * sa,
              deaths, be + ba, mu * (be + ba), se_in + sa_in]
        old_d, old_s = functional_components(y, s)
        old_raw = de + da
        y = [value + dt * flow for value, flow in zip(y, dy)]
        if min(y) < -1e-8 or not all(math.isfinite(v) for v in y):
            raise ArithmeticError('invalid state; reduce timestep or inspect parameters')
        new_d, new_s = functional_components(y, s)
        integral_d += dt * (old_d + new_d) / 2
        integral_s += dt * (old_s + new_s) / 2
        integral_raw_d += dt * (old_raw + y[2] + y[3]) / 2
        discrepancy = abs(sum(y[:6]) + y[6] - n0 - y[7])
        max_accounting_error = max(max_accounting_error, discrepancy)
        max_abnormal = max(max_abnormal, y[1] + y[3] + y[5])
        if retain_trajectory and (step % max(1, round(0.5 / dt)) == 0 or step == steps - 1):
            history.append({'time':(step + 1)*dt, 'state':y.copy(),
                            'functional_output':new_d + new_s})
    live = sum(y[:6])
    result = {
        'final_compartments':dict(zip(('P_E','P_A','D_E','D_A','S_E','S_A'), y[:6])),
        'live_cells':live, 'euploid_fraction':(y[0]+y[2]+y[4])/live,
        'abnormal_live_cells':y[1]+y[3]+y[5], 'max_abnormal_live_cells':max_abnormal,
        'cumulative_deaths':y[6], 'cumulative_divisions':y[7],
        'assumed_karyotype_changing_events':y[8], 'cumulative_senescence_transfers':y[9],
        'raw_differentiated_cell_time':integral_raw_d,
        'weighted_differentiated_cell_time':integral_d,
        'weighted_senescent_cell_time':integral_s,
        'weighted_functional_cell_time':integral_d+integral_s,
        'accounting_absolute_error':max_accounting_error,
    }
    if retain_trajectory:
        result['trajectory'] = history
    return result


def allele_sensitivity(alleles):
    products = alleles['allele_products']
    # Canonical consequence must be established by the actual sequence layer.
    def seq(key):
        p = products[key]
        return p.get('protein_sequence') or p.get('canonical_protein_sequence')
    if seq('L737Ter') != seq('L737Ter_N1002Lys_cis') or not seq('L737Ter'):
        raise ValueError('sequence layer does not establish masked cis missense')
    rows = []
    for phase, abundance, function, stop_abundance, stop_function in itertools.product(
            ('trans','cis'), (0.1,0.5,1.0), (0.0,0.5,1.0), (0.0,0.5), (0.0,0.5,1.0)):
        stop = stop_abundance * stop_function
        before = 0.5 * (abundance * function + stop) if phase == 'trans' else 0.5 * (1.0 + stop)
        after = 0.5 * (function + stop) if phase == 'trans' else before
        rows.append({'phase':phase,'missense_abundance':abundance,'missense_function':function,
                     'stop_abundance':stop_abundance,'stop_function':stop_function,
                     'functional_dose_before':before,'functional_dose_after_missense_abundance_restoration':after,
                     'change':after-before})
    return {'scope':'Algebraic sensitivity to unmeasured abundance/function. Additive model; dominant interference excluded. No dose-to-error calibration, NMD fraction estimate or drug response.',
            'rows':rows,
            'cis_changes':sorted(set(r['change'] for r in rows if r['phase']=='cis')),
            'zero_function_trans_changes':sorted(set(r['change'] for r in rows if r['phase']=='trans' and r['missense_function']==0)),
            'all_assumed_parameters':True}


def numerical_checks():
    scenarios = [Scenario(True,1.8,0.04,0.25,1.0,0.25),
                 Scenario(False,0.7,0.005,0.25,0.0,0.25),
                 Scenario(True,0.7,0.02,0.0,0.0,0.0,1.0)]
    worst = 0.0
    accounting = 0.0
    for s, duration in itertools.product(scenarios, HORIZONS):
        for arm in INTERVENTIONS:
            coarse = simulate(s,arm,dt=0.05,duration=duration)
            fine = simulate(s,arm,dt=0.025,duration=duration)
            for key in ('live_cells','raw_differentiated_cell_time','weighted_functional_cell_time'):
                worst = max(worst,abs(coarse[key]-fine[key])/max(1.0,abs(fine[key])))
            accounting = max(accounting,coarse['accounting_absolute_error'],fine['accounting_absolute_error'])
    if worst > 0.005 or accounting > 1e-6:
        raise ArithmeticError(f'convergence/accounting failed: {worst}, {accounting}')
    s = scenarios[0]
    assert simulate(s)==simulate(s,'functional_error_reduction',error_factor=1.0)
    s = replace(s,renewal=False)
    assert simulate(s)==simulate(s,'functional_error_reduction',error_factor=0.0)
    return {'max_relative_step_refinement_difference':worst,
            'max_absolute_accounting_error':accounting,
            'null_intervention_and_no_division_limits':'passed',
            'scope':'Numerical correctness only; not biological validation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alleles',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    start=time.monotonic()
    alleles=json.loads(args.alleles.read_text())
    checks=numerical_checks()
    rows=[]
    for duration, values in itertools.product(HORIZONS, itertools.product((True,False),(0.7,1.8),(0.005,0.02,0.04),
                                    (0.0,0.25),(0.0,1.0),(0.0,0.25))):
        s=Scenario(*values)
        baseline=simulate(s,duration=duration)
        for arm in INTERVENTIONS:
            r=baseline if arm=='untreated' else simulate(s,arm,duration=duration)
            r.update(intervention=arm,scenario=asdict(s),duration=duration,
                raw_output_ratio=r['raw_differentiated_cell_time']/baseline['raw_differentiated_cell_time'],
                weighted_output_ratio=r['weighted_functional_cell_time']/baseline['weighted_functional_cell_time'],
                euploid_fraction_change=r['euploid_fraction']-baseline['euploid_fraction'],
                required_differentiated_function_multiplier=max(0.0,
                    (baseline['weighted_functional_cell_time']-r['weighted_senescent_cell_time'])/r['weighted_differentiated_cell_time']))
            rows.append(r)
    illustrative={}
    for key,s in {'renewing':Scenario(True,1.8,0.04,0.25,1.0,0.25),
                  'nonrenewing':Scenario(False,0.7,0.02,0.25,0.0,0.25)}.items():
        illustrative[key]={'scenario':asdict(s),'arms':{a:simulate(s,a,retain_trajectory=True) for a in INTERVENTIONS}}
    # Alternative niche occupancy is a structural sensitivity, not measured biology.
    niche_sensitivity=[]
    for weight in (0.0,1.0):
        s=Scenario(True,0.7,0.02,0.25,0.0,0.25,weight)
        niche_sensitivity.append({'weight':weight,'baseline':simulate(s),
                                  'senolysis':simulate(s,'selective_senescent_removal')})
    payload={'scope':'Sequence-grounded allele consequences plus assumption-explicit population sensitivity. No whole-genome digital twin, fitted patient trajectory, drug efficacy, clinical probability, tumor model or first-founder endpoint.',
        'runtime':{'host':platform.node(),'machine':platform.machine(),'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start},
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'allele_input_sha256':hashlib.sha256(args.alleles.read_bytes()).hexdigest(),
        'time_unit':'arbitrary model-time unit, not patient days',
        'horizons':list(HORIZONS),'dt':0.05,'numerical_checks':checks,
        'competition_regimes':competition_regimes(),
        'parameter_status':'All population/intervention/function weights are assumptions. Error range overlaps Rozhok published SIMULATION inputs, not measured case rates. Initial 20% aneuploid fraction is not inferred from VAF or the retracted depth signal.',
        'allele_sensitivity':allele_sensitivity(alleles),
        'parameter_setting_count':96,'scenario_count':len(rows)//len(INTERVENTIONS),'rows':rows,
        'illustrative_trajectories':illustrative,'niche_sensitivity':niche_sensitivity,
        'limitations':['All abnormal karyotypes collapsed; no chromosome-specific fitness, allele-copy changes, LOH or malignant evolution.',
                       'Errors in already abnormal cells counted but do not change severity or create euploid reversions.',
                       'No direct N1002K or L737* function measurement; no fixed mapping from functional dose to error rate.',
                       'No fraction of grid points is a probability or a robustness claim outside this grid.',
                       'Finite assays may be inconclusive; a negative simulation does not refute a drug or mechanism in the child.']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'scenario_count':payload['scenario_count'],'runs':len(rows),
                      'elapsed_seconds':payload['runtime']['elapsed_seconds'],'checks':checks}))


if __name__=='__main__':
    main()
