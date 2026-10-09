#!/usr/bin/env node
/* Optional browser check; requires Playwright and its Chromium runtime. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const os = require('node:os');
const { pathToFileURL } = require('node:url');
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true, executablePath: process.env.EXPLORER_CHROMIUM });
  try {
    const page = await browser.newPage({ viewport: { width: 1360, height: 1100 } });
    const errors = [];
    const remoteRequests = [];
    page.on('pageerror', error => errors.push(error.message));
    page.on('request', request => { if (/^https?:/.test(request.url())) remoteRequests.push(request.url()); });
    await page.goto(pathToFileURL(path.join(__dirname, 'index.html')).href);
    assert.equal(await page.locator('#fraction-change').textContent(), '+16.8 pp');
    assert.equal(await page.locator('#output-change').textContent(), '−11.0%');
    assert.match(await page.locator('#output-values').textContent(), /51,111.0.*45,513.3/);
    const options = await page.evaluate(() => window.CELL_SELECTION_DATA.options);
    const arms = await page.locator('#arm option').evaluateAll(items => items.map(item => item.value));
    let checked = 0;
    for (const renewal of options.renewal) {
      await page.selectOption('#renewal', String(renewal));
      for (const division of options.aneuploid_division_ratio) {
        await page.selectOption('#division', String(division));
        for (const error of options.error_probability) {
          await page.selectOption('#error', String(error));
          for (const duration of options.duration) {
            await page.selectOption('#duration', String(duration));
            for (const arm of arms) {
              await page.selectOption('#arm', arm);
              const result = await page.evaluate(({ renewal, division, error, duration, arm }) => {
                const row = window.CELL_SELECTION_DATA.rows.find(row => row.duration === duration &&
                  row.scenario.renewal === renewal && row.scenario.aneuploid_division_ratio === division &&
                  row.scenario.error_probability === error && row.intervention === arm);
                return { row, fraction: document.getElementById('fraction-change').textContent,
                  output: document.getElementById('output-change').textContent };
              }, { renewal, division, error, duration, arm });
              const displayedNumber = text => Number(text.replace('−', '-').replace(/[^0-9.+-]/g, ''));
              assert.ok(Math.abs(displayedNumber(result.fraction) - result.row.euploid_fraction_change * 100) <= 0.050001);
              assert.ok(Math.abs(displayedNumber(result.output) - (result.row.raw_output_ratio - 1) * 100) <= 0.050001);
              checked++;
            }
          }
        }
      }
    }
    await page.locator('#reset').click();
    const downloadPromise = page.waitForEvent('download');
    await page.locator('#download').click();
    const download = await downloadPromise;
    const saved = JSON.parse(await fs.readFile(await download.path(), 'utf8'));
    const source = JSON.parse(await fs.readFile(path.join(__dirname, '..', 'population-results.json'), 'utf8'));
    const reference = source.illustrative_trajectories.renewing.arms.selective_aneuploid_killing;
    assert.equal(saved.selected.raw_differentiated_cell_time, reference.raw_differentiated_cell_time);
    assert.equal(saved.selected.euploid_fraction, reference.euploid_fraction);
    assert.equal(saved.selection.intervention, 'selective_aneuploid_killing');
    assert.match(saved.time_unit, /not patient days/);
    assert.equal(saved.source_sha256.length, 64);
    const screenshotDir = await fs.mkdtemp(path.join(os.tmpdir(), 'cell-selection-explorer-'));
    await page.screenshot({ path: path.join(screenshotDir, 'desktop.png'), fullPage: true });
    await page.setViewportSize({ width: 390, height: 844 });
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth), false);
    await page.screenshot({ path: path.join(screenshotDir, 'mobile.png'), fullPage: true });
    await page.locator('#duration').focus();
    await page.keyboard.press('8');
    await page.keyboard.press('Tab');
    assert.equal(await page.locator('#duration').inputValue(), '80');
    assert.match(await page.locator('#scenario-label').textContent(), /80 MODEL-TIME/);
    await page.locator('#reset').click();
    assert.equal(await page.locator('#fraction-change').textContent(), '+16.8 pp');
    assert.deepEqual(errors, []);
    assert.deepEqual(remoteRequests, []);
    console.log(JSON.stringify({ checked_comparisons: checked, default_example: 'passed',
      downloaded_source_values: 'exact match', keyboard_controls: 'passed',
      mobile_horizontal_overflow: false, browser_errors: errors,
      remote_requests: remoteRequests, screenshot_directory: screenshotDir }, null, 2));
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
