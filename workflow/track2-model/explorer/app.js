'use strict';
(() => {
  const data = window.CELL_SELECTION_DATA;
  const $ = id => document.getElementById(id);
  const number = value => value.toLocaleString('en-US', { maximumFractionDigits: 1, minimumFractionDigits: 1 });
  const signed = value => `${value > 0.00001 ? '+' : value < -0.00001 ? '−' : ''}${Math.abs(value).toFixed(1)}`;
  const descriptions = {
    untreated: 'The stored reference with no intervention. Both sides of the comparison are identical.',
    selective_aneuploid_killing: 'Adds 0.1 per model-time unit to the death rate of abnormal progenitors and differentiated cells. No chromosome correction is applied.',
    functional_error_reduction: 'Halves the assumed error probability during future divisions. Existing abnormal cells are not corrected.',
    stress_protection: 'Halves abnormal-cell death and senescence rates. It does not repair chromosome errors.',
    growth_slowing_only: 'Reduces progenitor division rates by 35%. All other rates stay at baseline.',
    growth_slowing_and_less_senescence: 'Reduces division rates by 35% and halves progenitor and differentiated-cell senescence rates.',
    selective_senescent_removal: 'Adds a senescent-cell removal rate of ln(2) during model intervals 5–6 and 15–16. These assumed pulses are not a clinical schedule.',
    senescent_removal_with_offtarget_loss: 'Uses the same senescent-cell pulses, plus an off-target loss hazard of −ln(1−0.02) for nonsenescent cells during each pulse.'
  };
  let selection;
  let baseline;
  let selected;

  function populate(id, options, label) {
    $(id).replaceChildren(...options.map(value => {
      const option = document.createElement('option');
      option.value = String(value);
      option.textContent = label(value);
      return option;
    }));
  }
  populate('division', data.options.aneuploid_division_ratio, n => `${n.toFixed(1)} × normal-cell rate`);
  populate('error', data.options.error_probability, n => `${number(n * 100)}% per division`);
  populate('duration', data.options.duration, n => `${n} model-time units`);
  $('source-hash').textContent = data.source_sha256;

  function barRow(label, value, max, pieces, sublabel = '') {
    const row = document.createElement('div');
    row.className = 'bar-row';
    const name = document.createElement('span');
    name.className = 'bar-label'; name.textContent = label;
    const track = document.createElement('div'); track.className = 'bar-track';
    track.setAttribute('aria-hidden', 'true');
    for (const [amount, className] of pieces) {
      const bar = document.createElement('span'); bar.className = `bar-fill ${className}`;
      bar.style.width = `${100 * amount / max}%`; track.append(bar);
    }
    const valueLabel = document.createElement('span'); valueLabel.className = 'bar-number';
    valueLabel.textContent = number(value);
    if (sublabel) { const sub = document.createElement('span'); sub.textContent = sublabel; valueLabel.append(sub); }
    row.append(name, track, valueLabel);
    return row;
  }

  function render() {
    selection = {
      renewal: $('renewal').value === 'true',
      aneuploid_division_ratio: Number($('division').value),
      error_probability: Number($('error').value),
      duration: Number($('duration').value),
      intervention: $('arm').value
    };
    const matching = data.rows.filter(row => row.duration === selection.duration &&
      ['renewal', 'aneuploid_division_ratio', 'error_probability'].every(key => row.scenario[key] === selection[key]));
    baseline = matching.find(row => row.intervention === 'untreated');
    selected = matching.find(row => row.intervention === selection.intervention);
    if (!baseline || !selected) throw new Error('This combination is absent from the retained grid.');
    const fractionChange = selected.euploid_fraction_change * 100;
    const outputChange = (selected.raw_output_ratio - 1) * 100;
    const isExample = Object.keys(data.default).every(key => selection[key] === data.default[key]);
    $('scenario-label').textContent = isExample ? 'RETAINED EXAMPLE' : `${selection.duration} MODEL-TIME UNITS`;
    $('arm-description').textContent = descriptions[selection.intervention];
    $('fraction-change').textContent = `${signed(fractionChange)} pp`;
    $('output-change').textContent = `${signed(outputChange)}%`;
    $('fraction-values').textContent = `${number(baseline.euploid_fraction * 100)}% untreated → ${number(selected.euploid_fraction * 100)}% selected`;
    $('output-values').textContent = `${number(baseline.raw_differentiated_cell_time)} → ${number(selected.raw_differentiated_cell_time)}`;
    let title = 'Compare composition with absolute output';
    let reading = 'A population fraction and accumulated output answer different questions. Read them together, alongside the absolute cell counts and the assumptions used here.';
    if (fractionChange > 0.01 && outputChange < -0.01) {
      title = 'A better fraction can accompany less output';
      reading = 'Here the chromosome-normal share rises while accumulated differentiated-cell output falls. The larger share alone would miss this cost. It does not demonstrate restoration of checkpoint function or correction of existing cells.';
    } else if (outputChange > 0.01) {
      title = 'More output under these particular assumptions';
      reading = 'The selected arm increases raw differentiated-cell time in this stored scenario. This is an assumed population response, not evidence of drug efficacy or measured muscle function. The chromosome-normal share remains a separate endpoint.';
    } else if (Math.abs(outputChange) <= 0.01) {
      title = 'Raw output is unchanged in this comparison';
      reading = selection.renewal ? 'Raw differentiated-cell time is unchanged. The optional weighted measure can still differ because it assumes contributions from senescent cells and impairment of neighboring cells.' : 'Without renewing progenitors, changing division or chromosome-error assumptions cannot create a benefit from future divisions. Other modeled cell losses can still affect output.';
    }
    if (selection.intervention === 'untreated') {
      title = 'The same reference on both sides';
      reading = 'Choose a hypothetical intervention to compare it with this untreated reference. Every comparison uses the same initial population and horizon.';
    }
    $('interpretation').textContent = title;
    $('reading').textContent = reading;
    const outputMax = Math.max(baseline.raw_differentiated_cell_time, selected.raw_differentiated_cell_time);
    $('output-bars').replaceChildren(
      barRow('Untreated', baseline.raw_differentiated_cell_time, outputMax, [[baseline.raw_differentiated_cell_time, 'reference']]),
      barRow('Selected', selected.raw_differentiated_cell_time, outputMax, [[selected.raw_differentiated_cell_time, '']])
    );
    const cellsMax = Math.max(baseline.live_cells, selected.live_cells);
    const counts = row => [row.live_cells - row.abnormal_live_cells, row.abnormal_live_cells];
    $('cell-bars').replaceChildren(...[baseline, selected].map((row, index) => {
      const [normal, abnormal] = counts(row);
      return barRow(index ? 'Selected' : 'Untreated', row.live_cells, cellsMax,
        [[normal, 'normal'], [abnormal, 'abnormal']], `${number(row.euploid_fraction * 100)}% normal`);
    }));
    const fields = [
      ['Chromosome-normal fraction (%)', row => row.euploid_fraction * 100],
      ['Living chromosome-normal cells', row => counts(row)[0]],
      ['Living abnormal cells', row => row.abnormal_live_cells],
      ['All living cells', row => row.live_cells],
      ['Raw differentiated-cell time', row => row.raw_differentiated_cell_time],
      ['Weighted functional cell-time (assumed)', row => row.weighted_functional_cell_time],
      ['Cumulative cell deaths', row => row.cumulative_deaths]
    ];
    $('result-table').replaceChildren(...fields.map(([label, get]) => {
      const tr = document.createElement('tr'); const th = document.createElement('th');
      th.scope = 'row'; th.textContent = label; tr.append(th);
      for (const row of [baseline, selected]) { const td = document.createElement('td'); td.textContent = number(get(row)); tr.append(td); }
      return tr;
    }));
  }

  function reset() {
    for (const [id, key] of [['renewal', 'renewal'], ['division', 'aneuploid_division_ratio'], ['error', 'error_probability'], ['duration', 'duration'], ['arm', 'intervention']]) {
      $(id).value = String(data.default[key]);
    }
    render();
  }
  for (const id of ['renewal', 'division', 'error', 'duration', 'arm']) $(id).addEventListener('change', render);
  $('reset').addEventListener('click', reset);
  $('download').addEventListener('click', () => {
    const payload = {
      schema_version: data.schema_version, title: data.title, method: data.method,
      source_file: data.source_file, source_sha256: data.source_sha256,
      model_source_sha256: data.model_source_sha256, scope: data.scope,
      time_unit: data.time_unit, parameter_status: data.parameter_status,
      selection, fixed_assumptions: data.fixed_assumptions,
      intervention_definition: descriptions[selection.intervention],
      interpretation: 'Aggregate compartments; not literal lineage tracking. Composition alone does not establish functional rescue.',
      comparison: {
        chromosome_normal_fraction_change_percentage_points: selected.euploid_fraction_change * 100,
        raw_differentiated_cell_time_change_percent: (selected.raw_output_ratio - 1) * 100
      },
      untreated: baseline, selected
    };
    const url = URL.createObjectURL(new Blob([JSON.stringify(payload, null, 2) + '\n'], { type: 'application/json' }));
    const anchor = document.createElement('a'); anchor.href = url;
    anchor.download = 'cell-selection-scenario.json'; anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  reset();
})();
