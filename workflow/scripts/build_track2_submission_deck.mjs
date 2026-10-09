import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const { SKILL_DIR, TMP_DIR, FINAL_PPTX, RUNTIME_PYTHON } = process.env;
for (const [name, value] of Object.entries({ SKILL_DIR, TMP_DIR, FINAL_PPTX, RUNTIME_PYTHON })) {
  if (!value || (name !== "RUNTIME_PYTHON" && !path.isAbsolute(value))) {
    throw new Error(`${name} must be set${name === "RUNTIME_PYTHON" ? "" : " to an absolute path"}`);
  }
}

const { resolvePresentationFont, finalizePresentation } = await import(
  pathToFileURL(path.join(SKILL_DIR, "container_tools/artifact_tool_utils.mjs")).href,
);

const W = 1280;
const H = 720;
const C = {
  paper: "#F7F4EC",
  white: "#FFFEFA",
  ink: "#12383D",
  text: "#203B3D",
  muted: "#596F70",
  teal: "#137F79",
  mint: "#76C7B8",
  paleMint: "#DCEDE6",
  gold: "#D69B42",
  paleGold: "#F3E6C9",
  coral: "#B85C48",
  paleCoral: "#F3DED6",
  rule: "#D1D8D2",
  dark: "#0E3339",
};
const FONT = resolvePresentationFont();
const presentation = Presentation.create({ slideSize: { width: W, height: H } });

function shape(slide, geometry, x, y, w, h, fill = "none", stroke = "none", width = 0, name) {
  return slide.shapes.add({
    geometry,
    ...(name ? { name } : {}),
    position: { left: x, top: y, width: w, height: h },
    fill,
    line: { fill: stroke, width },
  });
}

function rect(slide, x, y, w, h, fill, radius = 0, stroke = "none", width = 0, name) {
  const item = shape(slide, radius ? "roundRect" : "rect", x, y, w, h, fill, stroke, width, name);
  if (radius) item.borderRadius = radius;
  return item;
}

function text(slide, value, x, y, w, h, opts = {}) {
  const item = shape(slide, "textbox", x, y, w, h, "none", "none", 0, opts.name);
  item.text = value;
  item.text.style = {
    typeface: FONT,
    fontSize: opts.size ?? 24,
    bold: opts.bold ?? false,
    color: opts.color ?? C.text,
    alignment: opts.align ?? "left",
    verticalAlignment: opts.valign ?? "top",
    autoFit: opts.autoFit ?? "none",
    wrap: "square",
    insets: { top: 0, right: 0, bottom: 0, left: 0 },
    ...(opts.lineSpacing ? { lineSpacing: opts.lineSpacing } : {}),
  };
  return item;
}

function label(slide, value, x, y, w, color = C.teal) {
  return text(slide, value.toUpperCase(), x, y, w, 24, { size: 16, color, bold: true });
}

function pill(slide, value, x, y, w, fill, color = C.ink) {
  rect(slide, x, y, w, 32, fill, 15);
  text(slide, value, x + 12, y + 5, w - 24, 22, { size: 15, color, bold: true, align: "center" });
}

function frame(slide, title, n, dark = false) {
  slide.background.fill = dark ? C.dark : C.paper;
  const fg = dark ? C.white : C.ink;
  text(slide, title, 64, 43, 1060, 64, { size: 42, color: fg, bold: true });
  text(slide, `${String(n).padStart(2, "0")} / 08`, 1128, 53, 88, 26, {
    size: 16, color: dark ? C.mint : C.teal, bold: true, align: "right",
  });
  rect(slide, 64, 120, 86, 5, dark ? C.mint : C.teal);
  rect(slide, 64, 674, 1152, 1, dark ? "#36565A" : C.rule);
  text(slide, "jeweled  ·  MVA Track 2  ·  research proposal", 64, 686, 730, 20, {
    size: 14, color: dark ? "#BDD0CB" : C.muted,
  });
}

function note(slide, lines) {
  slide.speakerNotes.textFrame.setText(lines.join("\n"));
}

// 1 — Cover
{
  const slide = presentation.slides.add();
  slide.background.fill = C.paper;
  const coverPath = process.env.COVER_ART;
  if (!coverPath) throw new Error("COVER_ART must point to the extracted conceptual cover crop");
  const coverBytes = new Uint8Array(await fs.readFile(coverPath));
  slide.images.add({
    blob: coverBytes,
    contentType: "image/png",
    alt: "Conceptual abstract muscle-fiber illustration; not a patient image or experimental result",
    fit: "cover",
    position: { left: 625, top: 0, width: 655, height: H },
  });
  rect(slide, 0, 0, 640, H, C.paper);
  label(slide, "MVA  ·  Track 2", 64, 70, 300, C.teal);
  text(slide, "Preserving\nuseful muscle\nin MVA", 64, 126, 520, 216, { size: 55, color: C.ink, bold: true });
  text(slide, "Two distinct medicine hypotheses.\nOne chromosome-fidelity rule.", 66, 374, 495, 86, {
    size: 30, color: C.teal, bold: true,
  });
  text(slide, "Could existing medicines improve muscle function while preserving the cells needed for growth and repair?", 66, 482, 485, 76, {
    size: 22, color: C.text,
  });
  rect(slide, 64, 610, 490, 1, C.rule);
  text(slide, "EternaPeptix / jeweled  ·  laboratory proposal", 64, 628, 500, 24, { size: 16, color: C.muted, bold: true });
  text(slide, "No MVA drug rescue is established.", 64, 651, 490, 19, { size: 15, color: C.coral, bold: true });
  note(slide, [
    "Narration: This child has muscle wasting and mosaic variegated aneuploidy, a condition in which cells can carry different numbers of chromosomes. Our question is whether an existing medicine could improve muscle function while preserving the cells needed for growth and repair. We propose two medicines for laboratory testing.",
    "Sources: reports/jeweled_track2_report_submission.md; reports/appendices/track1-case-evidence.md.",
    "The cover crop reuses the prior conceptual AI-generated illustration. It is not a patient image, microscopy or experimental result.",
  ]);
}

// 2 — Candidate mechanism and uncertainty
{
  const slide = presentation.slides.add();
  frame(slide, "A BUBR1 supply hypothesis is testable", 2);
  text(slide, "A published MVA1 pattern supports testing usable BUBR1; it does not establish this child's allele function.", 64, 143, 1110, 42, {
    size: 22, color: C.muted,
  });

  rect(slide, 64, 211, 548, 245, C.white, 14, C.rule, 1);
  label(slide, "Observed candidate pair", 88, 232, 390);
  text(slide, "p.(Leu737Ter)", 88, 278, 250, 42, { size: 31, color: C.ink, bold: true });
  text(slide, "early stop; usable product may be reduced", 88, 322, 470, 32, { size: 20, color: C.text });
  rect(slide, 88, 372, 490, 1, C.rule);
  text(slide, "p.(Asn1002Lys)", 88, 391, 275, 36, { size: 28, color: C.teal, bold: true });
  text(slide, "stability and checkpoint function unmeasured", 88, 428, 470, 28, { size: 19, color: C.text });

  rect(slide, 636, 211, 580, 245, C.white, 14, C.rule, 1);
  label(slide, "Evidence boundary", 662, 232, 390);
  text(slide, "Selected MVA1 cells showed reduced BUBR1 abundance.", 662, 278, 510, 61, {
    size: 26, color: C.ink, bold: true,
  });
  text(slide, "Supports an abundance mechanism—not proof of drug response in this genotype.", 662, 354, 500, 58, {
    size: 21, color: C.text,
  });
  pill(slide, "PHASE STILL UNRESOLVED", 662, 423, 245, C.paleGold, C.ink);

  rect(slide, 64, 489, 1152, 112, C.paleMint, 14);
  label(slide, "Genetic stop rule", 90, 509, 250, C.teal);
  text(slide, "Trans is unobserved. If confirmed cis, retire this pair as the compound-heterozygous explanation.", 90, 543, 1085, 42, {
    size: 23, color: C.ink, bold: true,
  });
  note(slide, [
    "Narration: BUB1B makes the checkpoint protein BUBR1. Our stop-plus-missense pair matches a reported MVA1 pattern: reduced protein supply. The missense protein's stability and function remain unmeasured. Trans phase needs confirmation. Confirmed cis retires this pair as the compound-heterozygous explanation.",
    "Sources: reports/appendices/track1-case-evidence.md; reports/appendices/causal-interpretation.md; Suijkerbuijk et al. 2010, https://pubmed.ncbi.nlm.nih.gov/20516114/; Lischetti et al. 2020, https://pubmed.ncbi.nlm.nih.gov/33207204/.",
  ]);
}

// 3 — Niacin tests a proximal but incomplete restoration hypothesis
{
  const slide = presentation.slides.add();
  frame(slide, "Niacin tests usable BUBR1 restoration", 3);
  rect(slide, 64, 173, 548, 300, C.white, 14, C.rule, 1);
  rect(slide, 668, 173, 548, 300, C.white, 14, C.rule, 1);
  label(slide, "Human muscle engagement", 91, 201, 470, C.teal);
  text(slide, "Muscle NAD normalized in a small open-label mitochondrial-myopathy study.", 91, 245, 482, 72, { size: 27, color: C.ink, bold: true });
  text(slide, "Treated healthy controls also gained strength, so disease-specific benefit was not established.", 91, 343, 480, 66, { size: 21, color: C.text });
  label(slide, "Mechanism—two interventions apart", 695, 201, 480, C.teal);
  text(slide, "NMN increased BUBR1; genetic SIRT2 elevation improved selected outcomes in hypomorphic mice.", 695, 245, 480, 84, { size: 24, color: C.ink, bold: true });
  text(slide, "Neither result is niacin-mediated BUBR1 rescue.", 695, 354, 470, 48, { size: 21, color: C.text });
  rect(slide, 64, 505, 1152, 99, C.paleMint, 13);
  label(slide, "Untested bridge", 90, 524, 220, C.teal);
  text(slide, "Niacin → usable BUBR1 → functional muscle in this genotype.", 90, 557, 1080, 35, { size: 24, color: C.ink, bold: true });
  text(slide, "Prescription-niacin pediatric safety and effectiveness are unestablished.", 66, 625, 1144, 25, { size: 17, color: C.muted });
  note(slide, [
    "Narration: Niacin is our lead because it targets the protein-supply hypothesis. It raised muscle NAD in another myopathy. Separate experiments connect NAD and SIRT2 to BUBR1 stability. Niacin has not rescued this genotype. We test whether it increases usable protein and useful muscle.",
    "Niacin sources: Pirinen et al. 2020, https://pubmed.ncbi.nlm.nih.gov/32386566/; North et al. 2014, https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/; reports/appendices/track2-evidence-pipeline.md.",
    "Niacin does not inherit NMN or genetic SIRT2 effects. The current prescription-niacin label does not establish pediatric effectiveness.",
  ]);
}

// 4 — Apitegromab is a human-function comparator, not an MVA treatment
{
  const slide = presentation.slides.add();
  frame(slide, "Apitegromab tests muscle compensation", 4);
  text(slide, "Blocks promyostatin / latent-myostatin activation to support surviving muscle—not chromosome repair.", 64, 145, 1110, 45, { size: 22, color: C.muted });
  rect(slide, 64, 210, 548, 286, C.ink, 14);
  text(slide, "+1.8", 91, 237, 285, 102, { size: 90, color: C.mint, bold: true });
  text(slide, "HFMSE points vs placebo", 94, 340, 435, 31, { size: 24, color: C.white, bold: true });
  text(slide, "95% CI 0.30–3.32  ·  p=.019\nPooled-dose primary contrast", 94, 390, 454, 64, { size: 21, color: "#D8E6E0" });
  rect(slide, 668, 210, 548, 286, C.white, 14, C.rule, 1);
  label(slide, "SAPPHIRE · ages 2–12", 695, 237, 440, C.teal);
  text(slide, "Nonambulatory SMA on nusinersen or risdiplam.", 695, 274, 478, 51, { size: 23, color: C.ink, bold: true });
  rect(slide, 695, 340, 493, 69, C.paleGold, 9);
  text(slide, "20 mg/kg-only primary: p=.115 → gatekeeping stopped.", 709, 356, 464, 42, { size: 19, color: C.ink, bold: true });
  text(slide, "The 10 mg/kg estimate was exploratory. No MVA evidence.", 695, 428, 485, 44, { size: 18, color: C.text });
  rect(slide, 64, 530, 1152, 77, C.paleCoral, 13);
  text(slide, "US indication is SMA, not MVA. Fracture risk and recurring IV infusions matter.", 88, 553, 1090, 32, { size: 22, color: C.ink, bold: true, align: "center" });
  note(slide, [
    "Narration: Apitegromab is our muscle-support comparator. It blocks myostatin activation and improved motor function in a randomized spinal-muscular-atrophy trial. It now has an approved use, but no MVA evidence. Fractures are a labeled risk. We test useful muscle output, without claiming that larger muscles repair the chromosome defect.",
    "Sources: SAPPHIRE full report, https://orbi.uliege.be/bitstream/2268/339849/1/Safety%20and%20efficacy%20of%20apitegromab%20in%20nonambulatory%20type%202%20or%20type%203%20spinal%20muscular%20atrophy.pdf; protocol, https://cdn.clinicaltrials.gov/large-docs/20/NCT05156320/Prot_000.pdf; SAP, https://cdn.clinicaltrials.gov/large-docs/20/NCT05156320/SAP_001.pdf; US label, https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=b1cd8ab8-8a10-48ed-9806-75a44f22b963; reports/appendices/track2-candidate-shortlist.md.",
    "The pooled-dose primary contrast was positive; the prespecified 20 mg/kg-only primary comparison did not pass and stopped gatekeeping. Dose-specific/secondary estimates after that point were exploratory. This does not demonstrate MVA benefit.",
  ]);
}

// 5 — The key assay-design counterexample
{
  const slide = presentation.slides.add();
  frame(slide, "A better survivor fraction can still mean less tissue", 5, true);
  text(slide, "One simulated population scenario exposes the selection trap:", 64, 147, 1100, 34, { size: 23, color: "#D7E6DF" });
  text(slide, "+16.8", 84, 231, 470, 118, { size: 104, color: C.mint, bold: true });
  text(slide, "percentage points", 88, 351, 430, 34, { size: 24, color: C.white, bold: true });
  text(slide, "in chromosome-normal fraction", 88, 394, 465, 64, { size: 22, color: "#D7E6DF" });
  rect(slide, 626, 236, 4, 215, C.gold);
  text(slide, "−11.0%", 690, 231, 500, 118, { size: 100, color: "#F0B963", bold: true });
  text(slide, "total differentiated-cell time", 695, 351, 475, 38, { size: 24, color: C.white, bold: true });
  text(slide, "the model's measure of aggregate output", 695, 394, 475, 64, { size: 21, color: "#D7E6DF" });
  rect(slide, 84, 496, 1106, 82, "#19454B", 12);
  text(slide, "A favorable percentage can be produced by losing productive cells.", 108, 521, 1058, 38, { size: 26, color: C.white, bold: true, align: "center" });
  text(slide, "Assumption-based simulation—not patient data, a measured lineage result or a drug effect.", 86, 608, 1080, 25, { size: 17, color: "#C6D4CF" });
  note(slide, [
    "Narration: Our cell-selection explorer exposes a misleading result. In one simulated scenario, the chromosome-normal fraction rose sixteen point eight percentage points while total differentiated-cell time fell eleven percent. Other teams can inspect this tradeoff across declared assumptions. Our proposed experiment tracks starting cells, including those that die or stop dividing.",
    "Sources: workflow/track2-model/README.md; workflow/track2-model/population-results.json; workflow/track2-model/explorer/README.md; reports/appendices/founder-ledger.md.",
    "The displayed values are one assumption-based simulation scenario, not a patient measurement or drug response.",
  ]);
}

// 6 — Shared pass rule and candidate-specific target checks
{
  const slide = presentation.slides.add();
  frame(slide, "The experiment must show useful muscle", 6);
  text(slide, "Prespecify the same tissue-level rule for both mechanisms:", 64, 146, 1080, 34, { size: 23, color: C.muted });
  const x = [64, 458, 852];
  const fills = [C.paleMint, C.paleGold, C.paleCoral];
  const numbers = ["01", "02", "03"];
  const names = ["FUNCTION", "RENEWAL", "FIDELITY"];
  const copy = [
    "More force in established muscle OR more new-muscle output.",
    "Preserve the other functional cohort and renewable cell reserve.",
    "Chromosome errors do not worsen; track cells that stop or die.",
  ];
  for (let i = 0; i < 3; i += 1) {
    rect(slide, x[i], 209, 364, 243, C.white, 14, C.rule, 1);
    rect(slide, x[i], 209, 364, 11, i === 0 ? C.teal : (i === 1 ? C.gold : C.coral), 5);
    text(slide, numbers[i], x[i] + 24, 245, 82, 56, { size: 48, color: C.teal, bold: true });
    label(slide, names[i], x[i] + 24, 319, 302, i === 2 ? C.coral : C.ink);
    text(slide, copy[i], x[i] + 24, 352, 315, 79, { size: 21, color: C.text });
  }
  rect(slide, 64, 481, 1152, 111, C.ink, 13);
  label(slide, "Mechanism checks", 91, 500, 260, C.mint);
  text(slide, "Niacin: BUBR1 protein + function.   Apitegromab: target engagement in an active myostatin system.", 91, 534, 1080, 47, {
    size: 21, color: C.white, bold: true,
  });
  text(slide, "Tumor and cancer-treatment compatibility also needs review.", 90, 611, 1090, 26, { size: 17, color: C.muted });
  note(slide, [
    "Narration: We would compare affected and genetically corrected muscle models, using measured drug exposures. A medicine must improve either existing muscle or new muscle production, preserve the other, and avoid worsening chromosome errors. Independent confirmation must reproduce the result. The cancer history also requires checking that neither medicine helps tumor cells survive or interferes with cancer treatment.",
    "Sources: reports/appendices/track2-evidence-pipeline.md; reports/appendices/founder-ledger.md; reports/appendices/oncology-counterscreen.md.",
  ]);
}

// 7 — Search breadth, bounded conclusion
{
  const slide = presentation.slides.add();
  frame(slide, "A broad search supports a disciplined shortlist", 7);
  const metrics = [
    { x: 64, value: "15,357", label: "catalog/supplement entities", detail: "some supplement status unverified", color: C.teal },
    { x: 458, value: "40", label: "FDA product records", detail: "separate product-level layer", color: C.gold },
    { x: 852, value: "290", label: "route/compound appraisals", detail: "curated; not efficacy scores", color: C.coral },
  ];
  for (const m of metrics) {
    rect(slide, m.x, 164, 364, 194, C.white, 14, C.rule, 1);
    text(slide, m.value, m.x + 22, 181, 320, 66, { size: 54, color: m.color, bold: true });
    text(slide, m.label, m.x + 22, 257, 320, 31, { size: 19, color: C.ink, bold: true });
    text(slide, m.detail, m.x + 22, 296, 320, 43, { size: 16, color: C.muted });
  }
  rect(slide, 64, 382, 1152, 165, C.paleMint, 14);
  label(slide, "What the evidence leaves on the table", 90, 400, 560, C.teal);
  text(slide, "Niacin: restoration test.  NMN: mechanism benchmark.  Apitegromab: human-function comparator.", 90, 435, 1070, 50, {
    size: 22, color: C.ink, bold: true,
  });
  text(slide, "Everolimus stays conditional; common-compound follow-up changed no nominee and showed no BUBR1-dependent MVA muscle rescue.", 90, 489, 1080, 42, {
    size: 19, color: C.text,
  });
  text(slide, "Coverage is not exhaustive worldwide; retrieval and co-mentions are not efficacy assessments. Missing evidence stays unknown.", 66, 565, 1145, 43, {
    size: 17, color: C.muted,
  });
  note(slide, [
    "Narration: Our public drug census screened over fifteen thousand clinical and regulatory entries. It made apitegromab a stronger downstream candidate than everolimus for muscle function. Everolimus remains a conditional signaling probe. Other same-gene genetic and drug findings have explicit limitations. Missing evidence stays unknown rather than becoming a failed drug.",
    "Sources: workflow/track2-drug-census/README.md; workflow/track2-drug-census/decision-report.md; reports/appendices/track2-candidate-shortlist.md; reports/appendices/track2-route-decisions.md.",
    "The census retains 15,357 catalog/supplement entities and 40 FDA product entries in a separate layer; the ledger has 290 curated route/compound appraisals. These counts use different units and do not imply exhaustive market coverage or efficacy assessment of each entry.",
    "The common-compound extension at workflow/track2-drug-census/common-compound-extension.md is a bounded follow-up, not a complete supplement inventory; it did not establish BUBR1-dependent MVA rescue or change the provisional pair.",
  ]);
}

// 8 — The next work and a bounded scalability claim
{
  const slide = presentation.slides.add();
  frame(slide, "The next experiment depends on a qualified partner", 8, true);
  text(slide, "Model qualification", 84, 168, 298, 38, { size: 25, color: C.white, bold: true, align: "center" });
  text(slide, "→", 390, 166, 70, 43, { size: 36, color: C.gold, bold: true, align: "center" });
  text(slide, "Measured engagement", 464, 168, 318, 38, { size: 25, color: C.white, bold: true, align: "center" });
  text(slide, "→", 795, 166, 70, 43, { size: 36, color: C.gold, bold: true, align: "center" });
  text(slide, "Independent pass-rule test", 869, 168, 335, 38, { size: 24, color: C.white, bold: true, align: "center" });

  rect(slide, 84, 235, 1110, 133, "#19454B", 15);
  label(slide, "Needed next", 112, 260, 210, C.mint);
  text(slide, "Authorized myogenic model access + a laboratory partner.", 112, 299, 1020, 43, { size: 31, color: C.white, bold: true });
  text(slide, "No drug experiment has run. This is a research proposal—not treatment advice.", 86, 399, 1080, 32, { size: 20, color: "#D7E6DF", align: "center" });

  rect(slide, 84, 465, 1110, 105, "#19454B", 15);
  text(slide, "Reusable workflow", 111, 488, 245, 29, { size: 19, color: C.mint, bold: true });
  text(slide, "A public synthetic TRIP13 example did not inherit the BUBR1-restoration lead.", 111, 524, 1020, 34, { size: 22, color: C.white });
  text(slide, "Reuse demonstrates method portability—not a second-patient or treatment validation.", 86, 609, 1080, 27, { size: 17, color: "#C6D4CF", align: "center" });
  note(slide, [
    "Narration: We also ran the workflow on a public, synthetic TRIP13 example. It shows why a BUBR1 restoration hypothesis cannot simply transfer to another gene. That demonstrates reuse, not treatment validation. Next we need a consented muscle model and laboratory partner to test niacin and apitegromab separately.",
    "Sources: reports/jeweled_track2_report_submission.md; reports/appendices/track2-evidence-pipeline.md; workflow/track2-transfer/README.md.",
    "The public TRIP13 input is synthetic and the transfer example begins from selected coding variants, not raw-VCF discovery. Cell/lab access, staffing and costs remain unsecured.",
  ]);
}

await fs.mkdir(TMP_DIR, { recursive: true });
await fs.mkdir(path.dirname(FINAL_PPTX), { recursive: true });
const candidatePath = path.join(TMP_DIR, "candidate.pptx");
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);

const receiptPath = path.join(TMP_DIR, "jeweled_track2_slides_v4.validation.json");
const result = await finalizePresentation({
  workspaceDir: process.cwd(),
  candidatePath,
  finalPath: FINAL_PPTX,
  pythonExecutable: RUNTIME_PYTHON,
  integrityValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_package_integrity.py"),
  layoutValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_layout_geometry.py"),
  layoutArgs: ["--expected-slide-size-emu", "12192000,6858000", "--validate-bullet-geometry", "--validate-heading-fit"],
  explicitTotalSlideCount: 8,
  requiredNativeTableOwnerSlides: [],
  requiredNativeChartOwnerSlides: [],
  fontPolicy: { basis: "design", families: [FONT] },
  verifyArtifactToolImport: true,
  receiptPath,
});

for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i];
  const preview = await presentation.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(path.join(TMP_DIR, `slide-${i + 1}.png`), new Uint8Array(await preview.arrayBuffer()));
}
console.log(JSON.stringify({ finalPath: FINAL_PPTX, fontFamily: FONT, slides: presentation.slides.items.length, result }, null, 2));
