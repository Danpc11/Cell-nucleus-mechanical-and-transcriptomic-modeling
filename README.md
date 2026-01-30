# Cell nucleus mechanical and transcriptomic modeling

## Overview

This repository contains code and core results for an integrative analysis of liver fibrosis
progression in non-alcoholic fatty liver disease (NAFLD) and non-alcoholic steatohepatitis
(NASH). The study combines RNA-seq data from multiple human liver cohorts with
mechanically inspired models to link fibrosis progression, nuclear mechanics, and
transcriptional regulation in hepatocytes.

The central **hypothesis** is that progressive extracellular matrix (ECM) remodeling and tissue
stiffening during fibrosis lead to sustained nuclear deformation in hepatocytes, driving
stage-specific transcriptional programs that become dominant in advanced disease.

---

## Datasets

This project integrates three publicly available human liver RNA-seq datasets from the Gene
Expression Omnibus (GEO):

- [GSE130970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130970)
- [GSE135251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251)
- [GSE162694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162694)

All datasets are derived from snap-frozen human liver biopsies and span the full fibrosis
spectrum, including histologically normal controls (N) and fibrosis stages F1–F4. The
dominant disease etiology across cohorts is NAFLD/NASH-associated fibrosis.

---

### GSE130970 — NAFLD transcriptomic severity cohort

**Title:** Gene expression predicts histological severity and reveals distinct molecular
profiles of non-alcoholic fatty liver disease

- **Samples:** 78 human liver biopsies  
  - Normal controls: 6  
  - NAFLD/NASH: 72 (fibrosis stages F0–F4)
- **Species:** Homo sapiens  
- **Tissue:** Liver biopsy  
- **Sequencing:** Bulk RNA-seq  
- **Platform:** Illumina HiSeq  
- **Geographic origin:** United States
- **Age and sex:** Adult male and female patients; detailed per-sample age and sex metadata
  are reported in associated publications but are not uniformly encoded in the GEO series
  matrix.
- **Disease etiology:** NAFLD/NASH fibrosis

---

### GSE135251 — Multicenter NAFLD fibrosis cohort

**Title:** Transcriptomic profiling across the spectrum of non-alcoholic fatty liver disease

- **Samples:** 216 snap-frozen liver biopsies  
  - NAFLD/NASH: 206  
  - Healthy controls: 10
- **Species:** Homo sapiens  
- **Tissue:** Liver biopsy  
- **Sequencing:** Bulk RNA-seq  
- **Platform:** Illumina NextSeq 500  
- **Geographic origin:** Multicenter European  (United Kingdom)
- **Age and sex:** Adult male and female patients; age and sex distributions are available
  in the associated studies but are not consistently provided at the per-sample level in GEO.
- **Disease etiology:** Progressive NAFLD and NASH-associated fibrosis

---

### GSE162694 — Fibrosis-focused NASH cohort

**Title:** Molecular characterization and cell-type composition deconvolution of fibrosis in
NAFLD

- **Samples:** 143 human liver biopsies from NASH patients  
- **Species:** Homo sapiens  
- **Tissue:** Liver biopsy  
- **Sequencing:** Bulk mRNA-seq  
- **Platform:** Illumina-based RNA-seq  
- **Geographic origin:** United States
- **Age and sex:** Adult male and female patients; detailed demographic variables are
  reported in the original publication but are not fully structured in the GEO metadata.
- **Disease etiology:** NASH-associated fibrosis with emphasis on fibrosis severity

---

## Conceptual framework: hepatocyte nuclear mechanics across fatty liver fibrosis

### From fatty liver to nuclear mechanotransduction

This project is based on a hepatocyte-centered mechanobiological framework that links fatty
liver disease progression to nuclear deformation and transcriptional reprogramming.
 
#### **Healthy liver / Steatosis (N–F1)**

• Soft extracellular matrix 
• Low actomyosin tension 
• Rounded, compliant hepatocyte nucleus 
• Open chromatin organization 
• Predominantly metabolic transcriptional programs

#### **Progressive NASH (F2)**

• Initiation of ECM remodeling (collagen deposition, LOX activity) 
• Increased cytoskeletal tension and force transmission 
• Partial nuclear deformation 
• Chromatin reorganization and stress sensitivity 
• Mixed metabolic and inflammatory transcriptional programs 

#### **Advanced fibrosis (F3–F4)**

• Stiff, collagen-rich fibrotic matrix
• High mechanical load transmitted to the nucleus
• Flattened and strained hepatocyte nuclei
• Nuclear lamina stress and chromatin compaction
• Activation of mechanosensitive transcriptional programs
(ECM remodeling, cytoskeleton, YAP/TAZ–TEAD, survival pathways)


In early disease stages, hepatocytes experience lipid accumulation in a mechanically
compliant environment. Nuclear shape remains largely isotropic and gene expression is
dominated by metabolic programs.

As fibrosis progresses, ECM stiffening increases actomyosin-generated forces transmitted to
the nucleus via integrin–cytoskeleton–LINC complexes. At intermediate stages, chromatin
organization becomes mechanically sensitive.

## Mechanical modeling

To connect fibrosis progression with hepatocyte nuclear mechanics, we employ two
complementary modeling frameworks: (i) a weak power-law soft-matter rheological model and
(ii) generalized linear statistical models linking mechanical state to transcriptional
output.

---

### 1. Weak power-law soft-matter rheology

Liver tissue and hepatocyte nuclei are modeled as viscoelastic soft materials exhibiting
weak power-law behavior, a characteristic commonly observed in biological tissues.

The constitutive relation between stress and strain is expressed as:

σ(t) = G₀ · t^(−α) ⊗ ε(t), G(t) = G0 * t^(-alpha),   with 0 < alpha < 1

where:

- **σ(t)** is the mechanical stress as a function of time  
- **ε(t)** is the applied strain  
- **G₀** is an effective stiffness prefactor (generalized modulus)  
- **α** is the power-law exponent (0 < α < 1)  
- **⊗** denotes a temporal convolution  

This formulation implies that the material response is neither purely elastic (α = 0) nor
purely viscous (α = 1), but instead displays a continuous spectrum of relaxation times.

In the frequency domain, the complex shear modulus takes the form:

G*(ω) = G₀ · (iω)^α

with:

- **G*(ω)** complex shear modulus  
- **ω** angular frequency  
- **i** imaginary unit  

This representation captures the scale-free viscoelastic behavior of liver tissue and
hepatocyte nuclei during fibrosis progression.

---

### Biological interpretation

- **G₀ (effective stiffness)** increases with fibrosis stage due to extracellular matrix
  remodeling and collagen accumulation.
- **α (power-law exponent)** reflects the relative balance between elastic storage and
  viscous dissipation.
- Progressive fibrosis is associated with increasing **G₀** and reduced mechanical
  compliance, leading to enhanced force transmission to the nucleus.

The inferred mechanical state acts as a proxy for cumulative nuclear mechanical load.

---

### 2. Nuclear deformation and effective stress transmission

At the cellular level, actomyosin-generated forces are transmitted from the extracellular
matrix to the nucleus through integrins, the cytoskeleton, and the LINC complex. The
effective nuclear stress is modeled as:

**σ_nuc** = κ · **G_eff** · **ε_cell**

where:

- **σ_nuc** is the effective stress experienced by the nucleus  
- **κ** is a coupling coefficient capturing cytoskeleton–nucleus connectivity  
- **G_eff** is the effective tissue or cellular stiffness derived from the soft-matter model  
- **ε_cell** is the cellular strain  

This formulation links extracellular mechanical properties to nuclear deformation without
explicitly resolving intracellular force networks.

---

### 3. Generalized Linear Models (GLM) for transcriptomic coupling

To relate mechanical state variables to gene expression, we employ generalized linear
models of the form:

E[**y_g**] = **g⁻¹**(β₀ + β₁·S + β₂·M)

where:

- **y_g** is the expression level of gene *g*  
- **g⁻¹** is the inverse link function  
- **β₀** is the intercept  
- **S** is fibrosis stage (categorical or ordinal)  
- **M** is a mechanical proxy variable (e.g. G_eff or inferred σ_nuc)  
- **β₁, β₂** are regression coefficients  

For count-based RNA-seq data, the response is modeled using a negative binomial
distribution:

**y_g** ~ NB(**μ_g**, **θ_g**)

where:

- **μ_g** is the mean expression predicted by the linear predictor  
- **θ_g** is the dispersion parameter  

---

### Interpretation of GLM coefficients

- **β₁** captures stage-associated transcriptional effects independent of mechanics  
- **β₂** quantifies the contribution of mechanical load to gene expression  
- Genes with significant **β₂** values are interpreted as mechanically regulated or
  mechanosensitive  

This framework allows the identification of transcriptional programs that emerge
preferentially at advanced fibrosis stages (F3–F4), consistent with sustained nuclear
deformation.

---

### Summary

Together, the weak power-law rheological model and GLM-based transcriptomic coupling provide
a mechanistically interpretable bridge between fibrosis-driven tissue stiffening, nuclear
mechanics, and gene regulation in hepatocytes.
