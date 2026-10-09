#!/usr/bin/env python3
"""Finite abstract-context audit, with explicit manual exceptions; no efficacy ranking."""
import argparse,csv,hashlib,html,json,re
from collections import Counter,defaultdict
from pathlib import Path
from catalog import write_tsv,write_json,stamp
from literature_match import normalize
ROOT=Path(__file__).resolve().parent
# These annotations reflect inspection of the retained primary abstract, not full-paper
# verification. Exact compound roles and general context are intentionally separated.
MANUAL={
 '41700350':('lexical_false_match','Silver-Russell is a syndrome name in diagnostic testing. Silver was not administered.','none'),
 '39845481':('disease_model_intervention_other_lesion','Abstract reports NAC/ascorbic acid reduced segregation errors in Sod1-deficient oocytes. BUBR1 recruitment was reduced in this oxidative-stress model; actual-mutant MVA or muscle rescue was not tested.','NAC cross-lesion intervention; primary full text still needed'),
 '33798737':('disease_model_intervention_other_lesion','Abstract reports icariin improved motor coordination/learning in aged WT mice, with BUB1b among multiple upregulated markers and microbiota-transfer effects. BUB1B dependence and MVA rescue are not established.','icariin; unreviewed-primary lead'),
 '40273626':('disease_model_injury_exposure','Fentanyl and sevoflurane are exposures in the neonatal surgery/anesthesia injury model. Vitamin K2 is the reported protective intervention; BUB1B is one of multiple expression-network hubs. Exact K2 formulation and causal role need primary review.','vitamin K2; unreviewed-primary lead; do not promote matched anesthetic'),
 '33244299':('lexical_false_match','The tested compound is caffeic acid phenethyl ester (CAPE), not caffeic acid. Abstract reports neurodifferentiation/other-disease functional effects while BubR1 is attenuated. Parent-acid efficacy cannot be inherited from this ester.','CAPE; unreviewed-primary lead; derivative not equivalent'),
 '17464181':('lexical_false_match','The demethylating treatment is 5-aza-2-deoxycytidine (decitabine), not ordinary 2-deoxycytidine/doxecitine. Subphrase matching cannot transfer the result.','decitabine; covered by existing source appraisal'),
 '31611967':('lexical_false_match','Protein C-terminus binding is a gene-ontology term; it is not administration of human protein C concentrate.','none'),
 '25505175':('assay_reagent_or_background','Phenylalanine describes a sequence motif in BubR1, not a drug intervention.','none'),
 '14576056':('assay_reagent_or_background','Platelets are a measured hematologic phenotype of BubR1-mutant mice, not an administered platelet product.','none'),
 '37725942':('assay_reagent_or_background','Glucagon is an endogenous hormone measurement after genetic disruption of the insulin receptor/MAD2 interaction. The abstract does not report glucagon rescue.','none'),
 '32716490':('disease_model_injury_exposure','D-galactose was used as an aging/oxidative-stress challenge in Bub1b-heterozygous reproductive experiments, not a corrective medicine.','none'),
 '37151120':('genetic_rescue_with_drug_assay_reagent','Lifespan extension belongs to cleavage-resistant BubR1 flies. Colchicine is the checkpoint challenge used to measure kinetochore localization; the genetic result cannot be assigned to colchicine.','caspase-resistant BubR1 is genetic, not drug rescue'),
 '17272762':('assay_reagent_or_background','Acetylcholine is a vascular relaxation test stimulus in BubR1-mutant arteries; the response is impaired, not a systemic treatment rescue.','none'),
 '26883501':('assay_reagent_or_background','Fluorescein and Evans blue are permeability tracers in a BBB assay, including BubR1-hypomorphic mice. Their presence is not a therapeutic intervention.','none'),
 '24743044':('assay_reagent_or_background','Trypsin is used for protein digestion during proteomics. The pharmacological experiment uses tenovin-1 and reduces Bub-family proteins.','none'),
 '30250048':('cell_cycle_perturbation_or_harm','Lithium chloride and other GSK3 inhibitors shorten spindle-toxin arrest and impair checkpoint-protein recruitment/complex assembly. This is mechanistic perturbation, not a deficient-BUBR1 rescue result.','none'),
 '38604342':('cell_cycle_perturbation_or_harm','Simvastatin increases mitotic delay and abnormal cell division despite greater kinetochore BubR1. More checkpoint signal is not evidence of restored chromosome fidelity.','none'),
 '31053147':('disease_specific_clinical_care','Treosulfan is part of HSCT conditioning in an MVA1 patient with myelodysplasia. The abstract reports initial engraftment then graft rejection and fatal sepsis; no muscle/BUBR1 rescue.','none'),
 '31081598':('disease_specific_clinical_care','Dactinomycin was oncology treatment in BUB1B-associated bilateral Wilms tumor; severe adverse reaction stopped chemotherapy and nephrectomy followed. This is not correction of MVA muscle biology.','none'),
 '42330849':('cancer_cytotoxicity_or_resistance','Abstract reports griseofulvin binding Mad1/BubR1 and disrupting the checkpoint, with mitotic cancer-cell death. Norlichexanthone is a distinct TPR-targeting compound. Direct binding does not establish restoration.','griseofulvin direct-binding report; full primary assay review needed'),
 '38331723':('cross_disease_intervention_no_abstract_gene_link','Liraglutide improves viability/protein-degradation markers in high-glucose C2C12 cells. BUB1B is absent from the abstract; this does not establish a BUB1B-dependent myogenic rescue.','other-lesion muscle precedent only'),
 '35871686':('cross_disease_intervention_no_abstract_gene_link','Hesperetin is evaluated as a Cisd2-dependent intervention in naturally aged mice. The abstract does not establish BUB1B dependence or actual-mutant rescue.','CISD2 route, not direct BUB1B evidence'),
 '40999517':('cross_disease_intervention_no_abstract_gene_link','Fucoxanthin reduces apoptosis in oxidant-challenged placenta-derived MSCs. The abstract does not test BUB1B mutants or absolute muscle output.','oxidative-stress MSC precedent only'),
 '39009698':('cross_disease_intervention_no_abstract_gene_link','A niacinamide/hyaluronic-acid skin formula is studied in explants and adult skin. The abstract does not establish BUB1B restoration or muscle benefit.','skin formulation evidence does not transfer to niacin/MVA'),
 '39716891':('cross_disease_intervention_no_abstract_gene_link','CFI-402257/luvixasertib inhibits TTK in vascular injury/atherosclerosis models. The abstract does not establish deficient BUBR1 rescue or preserved myogenic progenitors.','other-lesion vascular intervention only'),
 '40197001':('lexical_false_match','The supplemented molecule is indole acetic acid, not acetic acid. AhR/Sirt2 fly longevity findings cannot be assigned to the parent subphrase.','indole acetic acid is a distinct compound'),
 '40167015':('cross_disease_intervention_no_abstract_gene_link','Navitoclax is a senolytic in aging-neurovascular experiments. Ganciclovir is used with a p16-3MR genetic sensitization system; it cannot inherit selective senolysis in unmodified MVA cells.','existing senolytic route; separate genetic trigger from pharmacologic drug'),
}
PER_DRUG={
 ('31787052','NIACINAMIDE'):('lexical_false_match','Nicotinamide occurs inside the name of NADPH oxidase; the abstract does not describe nicotinamide treatment.','none'),
 ('31787052','ANGIOTENSIN'):('disease_model_injury_exposure','Angiotensin II induces hypertension/fibrosis; low BubR1 attenuates this response. This does not show angiotensin as a treatment for BUBR1 deficiency.','none'),
 ('31787052','ANGIOTENSIN II'):('disease_model_injury_exposure','Angiotensin II is the hypertensive challenge in BubR1-low mice, not the rescuing intervention.','none'),
 ('24825348','NIACINAMIDE'):('lexical_false_match','The abstract names nicotinamide mononucleotide, a distinct precursor, and genetic SIRT2 overexpression. This lexical nicotinamide hit cannot itself establish niacinamide benefit; full-text niacinamide experiments are appraised separately.','NMN/SIRT2 existing appraisal'),
 ('19846757','NIACINAMIDE'):('lexical_false_match','The protective intervention is overexpression of nicotinamide phosphoribosyltransferase (NAMPT), not niacinamide administration.','genetic NAMPT, not niacinamide rescue'),
 ('19846757','NOCODAZOLE'):('disease_model_injury_exposure','Nocodazole induces tetraploidy to measure endothelial dysfunction; genetic NAMPT overexpression provides the protection.','none'),
 ('34027130','SIROLIMUS'):('lexical_false_match','Rapamycin occurs in the expansion of mechanistic target of rapamycin (mTOR). This observational oocyte transcriptome study states that no additional experimental work was performed.','none'),
 ('21369694','ETHYLENE'):('lexical_false_match','Ethylene describes a structural bridge in a combretastatin analog; ethylene is not a treatment.','none'),
 ('31771225','ACRIDINE'):('lexical_false_match','Acridine occurs within acridine-orange assay staining, not an acridine drug intervention.','none'),
 ('37274003','POLYETHYLENE GLYCOL'):('lexical_false_match','Polyethylene glycol diacrylate is a modified component of a surgical-guide resin, not ordinary PEG medication or MVA treatment.','none'),
 ('29044175','SIROLIMUS'):('lexical_false_match','Rapamycin appears inside mTORC1 terminology. The rescuing treatment in Pompe myocytes is recombinant GAA enzyme, not sirolimus.','none'),
}
# Additional bounded abstract checks for plausible checkpoint or clinical leads.
MANUAL.update({
 '24695965':('assay_reagent_or_background','Rapamycin is a heterodimerization trigger for engineered FKBP-MAD1/FRB-MIS12 recruitment. Restored checkpoint activity depends on the artificial constructs; it is not demonstrated drug rescue in unmodified MVA cells.','none'),
 '39737806':('cell_cycle_perturbation_or_harm','GSK923295 inhibits CENP-E in spermatogonial experiments. CENP-E disruption is associated with misalignment, impaired self-renewal and aneuploidy; this is not a BUBR1-restoring intervention.','none'),
 '38856693':('cancer_cytotoxicity_or_resistance','Prexasertib targets SF3B1-mutant hematopoietic stem cells with BUBR1/CDC27 mis-splicing and abrogates engraftment. Selective elimination of neoplastic cells is not restoration of MVA progenitors.','none'),
 '41579174':('disease_model_injury_exposure','Gentamicin induces nephrotoxicity in the analyzed rat datasets. BUB1B is part of a response-network core, not evidence that gentamicin corrects BUB1B deficiency.','none'),
 '42667124':('assay_reagent_or_background','Chorionic gonadotropin is described as a secreted trophoblast hormone. Kinase perturbations study placental cell differentiation; no hCG-mediated BUB1B rescue is reported.','none'),
 'PPR1286781':('assay_reagent_or_background','Preprint of trophoblast-kinase work: chorionic gonadotropin is a secreted hormone marker, not an established BUB1B-rescue treatment.','none'),
 '42278672':('computational_or_biomarker_association','BUB1B is among exploratory blood markers associated with dupilumab response in atopic dermatitis. The study does not establish deficient BUB1B restoration or MVA benefit.','none'),
 '21550059':('disease_model_injury_exposure','Angiotensin II/oxidant exposure induces BubR1 expression in young vascular smooth-muscle cells. This is a signaling-stimulus experiment, not demonstrated mutant skeletal-muscle rescue.','none'),
})
PER_DRUG.update({
 ('22886373','CHROMIUM'):('lexical_false_match','The exposure is toxic hexavalent chromium, not evidence supporting a medicinal chromium formulation. Oxidation state and compound identity cannot be ignored.','none'),
 ('22886373','ACETYLCYSTEINE'):('disease_model_intervention_other_lesion','NAC reduces chromium-induced ATP depletion and cell-cycle arrest in hepatocytes. The checkpoint-protein rescue described for a p53 inhibitor cannot automatically be assigned to NAC. No BUB1B-mutant muscle experiment.','existing NAC route; toxicant-stress precedent only'),
})
# Final finite pass: all 28 formerly unresolved gene-literal rows had abstracts.
# These are manual abstract-role dispositions, not full-text efficacy adjudication.
FINAL_PASS_RECORDS={
 'IND608907204','39083874','35217622','39687194','34919076','36389792',
 '35326770','IND606941299','32078082','25567764','18723675','20646066',
 '20372054','KR20070003777','CN1871349','US2007238098','WO2005045029',
 'EP1688492','11864974','15870702','15220334',
}
MANUAL.update({
 'IND608907204':('disease_model_injury_exposure','Busulfan conditions immunodeficient mice for human-PBMC engraftment and a graft-versus-host liver-injury model. BUB1B is part of an activated T-cell expression signature, not a corrected BUB1B lesion. Same study also indexed as PMID39687194.','none'),
 '39687194':('disease_model_injury_exposure','Busulfan is conditioning for a xenogeneic graft-versus-host model, not a BUB1B-restoring treatment. Liver injury and activated human T-cell expression are measured. AGR record IND608907204 duplicates this study.','none'),
 '39083874':('disease_model_injury_exposure','Propylparaben exposure disrupts mouse oocyte maturation, spindle attachment and chromosome alignment, with persistent BubR1/BUB3 checkpoint activation. The experiment reports harm, not correction of an inherited BUB1B defect.','none'),
 '35217622':('assay_reagent_or_background','ATP supplies bond hydrolysis in mitotic extracts/reconstituted ubiquitylation reactions for UBR5-dependent checkpoint-complex disassembly. This is biochemical substrate use, not administration of an ATP medicine.','none'),
 '34919076':('computational_or_biomarker_association','Dexamethasone, nicotine and resveratrol are nominated by gene-chemical network analysis for disc degeneration. The abstract reports no drug-treatment rescue experiment or BUB1B-mutant intervention.','none'),
 '36389792':('computational_or_biomarker_association','Drug names are ranked candidates from blood-transcriptome intersections of HIV and COVID-19 datasets. The abstract reports no testing of those compounds, BUB1B restoration or MVA efficacy.','none'),
 '35326770':('disease_model_injury_exposure','Polymyxin B is tested for lung epithelial toxicity, including cell death and an altered gene network containing BUB1B. Network membership does not make this a BUB1B-rescue intervention.','none'),
 'IND606941299':('assay_reagent_or_background','Putrescine is required culture-medium supplementation for arginase-deficient CHO cells. Polyamine starvation arrests growth and alters checkpoint-gene transcripts; no BUB1B-mutant rescue was tested. Same study also indexed as PMID32078082.','none'),
 '32078082':('assay_reagent_or_background','The experiment deprives arginase-deficient CHO cells of required polyamines and measures growth/transcript changes. Putrescine is a culture requirement, not a demonstrated MVA therapy. AGR record IND606941299 duplicates this study.','none'),
 '25567764':('assay_reagent_or_background','The tested agent is amyloid-beta peptide. Vinblastine identifies a tubulin-binding site near the peptide interaction, not a reported vinblastine rescue treatment.','none'),
 '18723675':('cell_cycle_perturbation_or_harm','Bleomycin induces mitotic DNA damage; UCN-01 inhibits CHK1 and reduces mitotic exit. RPA phosphorylation/genetic states determine the response. This is damaged-mitosis biology, not drug rescue of BUB1B-deficient muscle.','none'),
 '20646066':('cancer_cytotoxicity_or_resistance','Curcumin perturbs microtubules, accumulates checkpoint proteins including BubR1 at kinetochores, delays mitosis and induces MCF-7 apoptosis. Vinblastine/paclitaxel are combination comparators. Checkpoint activation here accompanies cancer-cell inhibition, not restoration.','none'),
 '20372054':('genetic_rescue_with_drug_assay_reagent','Actinomycin D displaces the chromosomal passenger complex and causes misalignment, cytokinesis failure and checkpoint override. Engineered Cenp-B-INCENP recruitment restores selected checkpoint outputs but not alignment/cytokinesis; that rescue cannot be assigned to the drug.','none'),
 '15870702':('cancer_cytotoxicity_or_resistance','Doxorubicin induces hepatoma-cell death through apoptosis or mitotic catastrophe depending on exposure. BubR1 loss and phosphatidylserine externalization are measured effects, not evidence of MVA correction.','none'),
 '15220334':('cell_cycle_perturbation_or_harm','Distamycin A (stallimycin) alters scaffold/matrix and centromeric domains, increases BubR1, and causes mitotic chromosomal missegregation. More BubR1 in this context is not evidence of restored fidelity or mutant rescue.','none'),
})
for patent_id in ['KR20070003777','CN1871349','US2007238098','WO2005045029','EP1688492']:
    MANUAL[patent_id]=('patent_risk_prediction_unvalidated','Patent abstract proposes CYP2C8/BUB1B SNP testing to estimate paclitaxel-associated granulocytopenia. It is a diagnostic-risk claim, not a BUB1B-rescue experiment; the abstract does not establish clinical validation. Related patent publications must not count as independent experiments.','clinical validity UNKNOWN from abstract')
PER_DRUG.update({
 ('11864974','OPIUM'):('lexical_false_match','Opium denotes the natural source class of purified noscapine. The experiment tests noscapine, not opium; source material and isolated alkaloid are not equivalent drugs.','none'),
 ('11864974','NOSCAPINE'):('cell_cycle_perturbation_or_harm','Noscapine changes microtubule dynamics, removes kinetochore tension and arrests HeLa cells in mitosis. Checkpoint localization changes are studied; the abstract does not show rescue of BUB1B-mutant tissue.','none'),
 ('15870702','PHOSPHATIDYL SERINE'):('assay_reagent_or_background','Phosphatidylserine externalization is measured as a membrane/cell-death marker after doxorubicin exposure. It is not an administered phosphatidylserine treatment.','none'),
})
CANCER=re.compile(r'\b(cancer|carcinoma|tumou?r|leukemia|leukaemia|lymphoma|glioblastoma|melanoma|neuroblastoma|myeloma|sarcoma)\b',re.I)
COMPUTE=re.compile(r'\b(docking|connectivity map|bioinformatics|bioinformatic|in silico|virtual screening|drug repositioning|drug repurposing)\b',re.I)
WET=re.compile(r'\b(treated|treatment with|administered|xenograft|knockdown|knockout|western blot|cell viability|cytotoxic|in vitro|in vivo)\b',re.I)
CELL=re.compile(r'\b(mitotic|spindle|microtubule|kinetochore|cell cycle|meiotic|oocyte|chromosome segregation|chromosomal)\b',re.I)
BENEFIT=re.compile(r'improv|ameliorat|alleviat|rescu|attenuat|restor|protect|revers|extend',re.I)
RESTORATIVE=re.compile(r'life.?span|longevi|senescen|aging|ageing|muscle|progeri|sarcopen|frailty|neurogenesis|cognitive|myogen|regener',re.I)

def plain(x):return ' '.join(html.unescape(re.sub('<[^>]+>',' ',x or '')).split())
def context(row,abstract):
    rid=row['record_id'];drug=row['drug_name'];nk=normalize(abstract);aliases=json.loads(row['matched_normalized_aliases'])
    if (rid,drug) in PER_DRUG:return *PER_DRUG[rid,drug],'manual_abstract_review'
    if rid in MANUAL:return *MANUAL[rid],'manual_abstract_review'
    if 'cofactor' in aliases and drug=='FOLITIXORIN':
        return 'lexical_false_match','Cofactor is an ordinary biological term and an ambiguous research-code synonym; no folitixorin intervention is established.','none','manual_lexical_review'
    if not abstract:return 'unresolved_abstract_missing','No abstract is available; the title alone is insufficient to classify the compound role.','abstract/full text needed','explicit_missing_data'
    # A lower-priority, gene-only-in-full-text record is not elevated from generic beneficial language.
    if row['direct_gene_name_in_title_or_abstract']!='True':
        return 'unresolved_gene_context_requires_full_text','Drug alias appears in title/abstract, but the gene literal does not. Establish the gene connection and intervention identity from primary text.','unresolved','scope_rule'
    # Context tags below are transparent routing rules, not drug efficacy judgments.
    if COMPUTE.search(abstract) and not WET.search(abstract):
        return 'computational_or_biomarker_association','Abstract reports computational/association methods without an identified treatment experiment under this routing rule. No efficacy inference.','none','abstract_context_rule'
    if CANCER.search(abstract):
        return 'cancer_cytotoxicity_or_resistance','Abstract has an oncology context. The matched compound may be treatment, comparator, resistance exposure or prediction; exact role and direction need primary review.','none','abstract_context_rule'
    if drug.lower() in {'eosin','trypan blue free acid','fluorescein','evans blue free acid','trypsin, crystallized'}:
        return 'assay_reagent_or_background','Matched compound is a conventional staining/tracer/proteomics reagent in this abstract context; a therapeutic role is not established.','none','abstract_context_rule'
    if CELL.search(abstract) and drug in {'NOCODAZOLE','PACLITAXEL','COLCHICINE','VINBLASTINE','DOCETAXEL ANHYDROUS','DOXYCYCLINE ANHYDROUS','TETRACYCLINE','DOXORUBICIN','ETOPOSIDE','TRICHOSTATIN'}:
        return 'assay_reagent_or_cell_cycle_perturbation','Abstract describes cell-cycle/checkpoint research with a drug or induction-system alias. This does not by itself show restoration of deficient BUBR1 or useful muscle.','none','abstract_context_rule'
    if COMPUTE.search(abstract):
        return 'computational_or_biomarker_association','Computational/association methods appear in the abstract, potentially with other experiments; compound-specific validation needs full text.','none','abstract_context_rule'
    return 'unresolved_needs_full_text','Available abstract does not securely establish the matched compound role or its relevance to BUB1B restoration.','unresolved','explicit_unresolved'

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--corpus',required=True,type=Path);ap.add_argument('--out',default=ROOT,type=Path);a=ap.parse_args()
    raw=a.corpus.read_bytes();data=json.loads(raw);ci={(r['source'],r['id']):r for r in data['results']}
    mentions=list(csv.DictReader((a.out/'literature-mentions.tsv').open(),delimiter='\t'));rows=[];candidates=[]
    for row in mentions:
        r=ci[row['source'],row['record_id']];abstract=plain(r.get('abstractText',''))
        cat,note,lead,method=context(row,abstract)
        out={k:row[k] for k in ['catalog_id','drug_name','source','record_id','pmid','doi','title','matched_normalized_aliases','source_url']}
        out.update({'gene_literal_in_title_abstract':row['direct_gene_name_in_title_or_abstract'],'abstract_available':bool(abstract),'triage_category':cat,'triage_method':method,'abstract_supported_note':note,'follow_up':lead,'full_primary_verification':'NOT_COMPLETED_BY_THIS_ABSTRACT_AUDIT','actual_BUB1B_variant_muscle_rescue':'NOT_ESTABLISHED_BY_THIS_ROW','abstract_sha256':hashlib.sha256(abstract.encode()).hexdigest()})
        rows.append(out)
        if (RESTORATIVE.search(abstract) and BENEFIT.search(abstract) and not CANCER.search(plain(r.get('title','')))) or 'unreviewed-primary lead' in lead:
            candidates.append(out)
    write_tsv(a.out/'literature-triage.tsv',rows,list(rows[0]));write_tsv(a.out/'literature-restorative-context-audit.tsv',candidates,list(rows[0]))
    priority=[r for r in rows if r['gene_literal_in_title_abstract']=='True']
    summary={'built_utc':stamp(),'corpus_sha256':hashlib.sha256(raw).hexdigest(),'co_mention_rows':len(rows),'gene_literal_priority_rows':len(priority),'gene_literal_priority_unique_records':len({(r['source'],r['record_id']) for r in priority}),'all_categories':dict(Counter(r['triage_category'] for r in rows)),'priority_categories':dict(Counter(r['triage_category'] for r in priority)),'priority_methods':dict(Counter(r['triage_method'] for r in priority)),'restorative_context_rows':len(candidates),'restorative_context_unique_records':len({(r['source'],r['record_id']) for r in candidates}),'scope':'All specific alias co-mentions receive transparent abstract-context routing. Targeted manual inspection covers restorative/same-gene/injury/identity-confusion cases. This is a finite completion audit, not full-text adjudication of all papers. Context-rule categories are not article-level efficacy conclusions; unresolved stays unresolved.','additional_primary_followups':['icariin PMID33798737: aged WT mice; BUB1b is one of many changed markers','vitamin K2 PMID40273626: neonatal anesthetic injury; exact formulation and BUB1B mediation unresolved','CAPE PMID33244299: distinct ester, not caffeic acid; BubR1 attenuated in differentiation','NAC PMID39845481: Sod1-deficient oocyte segregation, existing route/cross-lesion evidence','griseofulvin PMID42330849: reported Mad1/BubR1 binding and cancer cytotoxicity; no restoration claim'],'no_new_efficacy_ranking':True,'final_unresolved_pass':{'originally_unresolved_rows':28,'source_records':21,'available_abstracts':21,'missing_abstracts':0,'status':'All 28 received manual abstract-role dispositions; no new restorative intervention identified in this final set. Full-primary clinical validity and mutant efficacy remain unestablished.','record_ids':sorted(FINAL_PASS_RECORDS)}}
    write_json(a.out/'literature-triage-summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
