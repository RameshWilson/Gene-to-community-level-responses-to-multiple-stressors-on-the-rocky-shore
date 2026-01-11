# Brighton warming × sewage-pollution experiment

Code and data for the manuscript: **Gene-to-community level responses to multiple stressors on the rocky shore**.

## Repository structure
- `Brighton project.Rproj`: RStudio project file (sets the working context).
- `Key analyses/`: core manuscript analysis scripts.
- `Additional analyses/`: diagnostics, sensitivity checks, and supplementary figure/processing scripts.
- `Dataframes/`: CSV inputs used by the scripts and selected large supplementary outputs (see below).

## Quick start (recommended workflow)
1. Install **Git**.
2. Clone this repository.
3. Open `Brighton project.Rproj` in RStudio.
4. Install R packages as needed (each script lists required packages at the top).
5. Run scripts from within the project.

> Note: scripts use the `{here}` package for portable paths (e.g. `here("Dataframes", "Full Brighton dataframe.csv")`).  
> If needed: `install.packages("here")`.

## Requirements
- R (analyses run on **R 4.4.3**).
- RStudio (recommended).
- Git.

---

## Transcriptomics

### Raw RNA-seq data
- Raw sequencing data are deposited at **ArrayExpress** (Accession: **E-MTAB-16541**).

### Derived transcriptomics outputs included in this repository 
To support transparency and interpretation, two large workbooks are included in the 'Dataframes' folder:
- **Annotated DEG list and interaction classes**: DEG-level outputs across contrasts; log2FC, significance, annotations/database hits, and interaction classifications 
- **GO enrichment results**: GO term enrichment summaries for DEG sets across contrasts

### Illustrative DESeq2 contrast snippet
- Provided in 'Additional analyses' folder
- This is a **minimal script** showing the DESeq2 model design and contrast definitions used for Figure 4 of the manuscript/summary tables above
- It is **not** a full transcriptomics pipeline: see ArrayExpress for raw data provenance and workflow context.

---

## SI figure: precipitation × sewage releases

This project includes code to generate a **supplementary figure** that visually compares:
- daily precipitation (Peacehaven area), and
- daily sewage release durations (Newhaven outfall).

### Sewage release data (Southern Water Rivers & Seas Watch)
Sewage release event history was downloaded from Southern Water’s **Rivers and Seas Watch** portal (publicly accessible; downloadable table).
> Note: even though the portal is publicly accessible, the dataset is third-party governed; see `DATA_LICENSE.md` for attribution/reuse guidance.

### Precipitation data (Visual Crossing — not redistributed here)
Daily precipitation was obtained from **Visual Crossing** during a time-limited paid subscription.

**Public repository policy**
- This repository **does not redistribute** Visual Crossing raw downloads (e.g. daily CSV time series), in line with Visual Crossing terms and conditions: https://www.visualcrossing.com/weather-services-terms/
- The SI figure script is included **for transparency**, but it **cannot run from this repository alone** unless one supplies their own precipitation file.

**To fully reproduce the SI figure**
1. Obtain your own daily precipitation time series for **01/06/2023–15/09/2023** for Peacehaven (or equivalent location) via Visual Crossing’s query builder (or an alternative provider).
2. Save the downloaded file locally (not committed) as:
   - `Dataframes/Precipitation data.csv`
3. Run:
   - `Additional analyses/Precipitation x sewage releases_Additional.R` (noting variable naming conventions within the script)

---

## Reproducibility notes
- Paths are project-relative via `{here}`.
- Figures are produced interactively unless a script explicitly writes outputs.

## How to cite
If you use this repository, please cite:
- **Manuscript:** citation + DOI (to be added)
- **Data**: https://doi.org/10.5281/zenodo.18121768

## Licence
- Code licence: see `LICENSE` (MIT).
- Repository data licence: see `DATA_LICENSE.md` (CC BY 4.0 for author-generated project datasets and derived tables), **excluding** third-party governed inputs noted there.
- Third-party data (Visual Crossing; Southern Water Rivers & Seas Watch) are governed by their respective terms (not covered by this repo’s licences).

## Contact
- Lead author: Ramesh Wilson  
- Email: rameshwilson14@gmail.com  
- GitHub: `RameshWilson`
