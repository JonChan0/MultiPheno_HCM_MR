# Prioritising modifiable risk factors and biomarkers in hypertrophic cardiomyopathy: a Mendelian randomization approach

**Authors:** Jonathan H Chana<sup>a,b</sup>, Christopher Grace<sup>a,b</sup>, Mohsen Mazidi<sup>c</sup>, Robert Clarke<sup>c</sup>, Hugh Watkins<sup>a,b</sup>, Anuj Goel<sup>a,b</sup>

**Affiliations:**
* <sup>a</sup> Division of Cardiovascular Medicine, Radcliffe Department of Medicine, University of Oxford, UK
* <sup>b</sup> Centre for Human Genetics, Nuffield Department of Medicine, University of Oxford, UK
* <sup>c</sup> Nuffield Department of Population Health, University of Oxford, UK

---

## Introduction

Hypertrophic cardiomyopathy (HCM) is one of the most common inherited heart diseases. It is responsible for driving heart failure and sudden cardiac death in patients. Numerous phenotypes have been associated with HCM, including:
* Clinical risk factors such as BMI
* Plasma proteins such as N-terminal pro-brain natriuretic peptide (NTproBNP)
* Cardiac abnormalities such as fibrosis

These associated phenotypes allow for the predictive modelling of incident HCM disease, but constructing such a model requires appropriate feature selection. As such, this study evaluated
such phenotypes for causal association with HCM using Mendelian randomization (MR). The objective was to assess whether a causality-based prioritization of upstream risk factors improved predictive modelling.

---

## Methods

* Bidirectional, two-sample MR was performed using the inverse-variance-weighted method or the Wald ratio estimate (for multiple and single instrument(s), respectively).
* Instruments were selected from lead variants for genome-wide significant loci for each phenotype.
* Sensitivity analyses, including MR-PRESSO, MR-RAPS, and weighted median, were applied to validate signal robustness.
* Significant causal associations were investigated by colocalisation with HCM via pairwise `coloc` to clarify the loci underpinning these causal relationships.
* Cox regression modelling with age as the timescale was applied to 548 incident cases and 466,078 controls from the UK Biobank. This was done after the left-truncation of prevalent cases and right-censoring.
* Model performance was evaluated by validation-fold C-index and time-dependent AUC. This used a 5-fold cross-validation scheme with bootstrapping (n=100) and was benchmarked against LASSO and Elastic Net Cox regression.

---

## Results

* MR analyses identified numerous significant associations between HCM-linked phenotypes and the disease.
* Causal associations with HCM were replicated for:
    * BMI (OR per SD = 1.60 [1.34 – 1.90] 95% CI)
    * Diastolic blood pressure (OR per SD = 1.05 [1.02 – 1.07] 95% CI)
* These findings prioritize hypertension and obesity as upstream risk factors.
* Novel causality was also identified between HCM and plasma protein biomarkers:
    * NTproBNP (SD per doubling of odds(HCM) = 0.10 [0.06-0.15] 95% CI)
    * Troponin T (SD per doubling of odds(HCM) = 0.06 [0.02-0.10] 95% CI)
* Colocalisation analyses further implicated `MYOZ1` and `BAG3` loci as those containing causal variants that underpin these associations.
* Cox regression modelling demonstrated no significant performance gain upon using MR as a basis for feature prioritisation.

---

## Conclusions

* This study applied MR to clarify causal relationships between HCM and both its potential upstream drivers and its downstream biomarkers.
* This process aided in the prioritization of modifiable risk factors like BMI and identified causality with plasma protein biomarkers such as NTproBNP.
* Such genetic evidence, however, does not translate to stronger predictive utility for incident HCM disease.

![ASHG Poster](https://github.com/JonChan0/MultiPheno_HCM_MR/blob/main/docs/ASHG_Poster.png)