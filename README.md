# Causal Role of IL-33/ST2 in Dementia: A Mendelian Randomisation Study

### Overview
This project investigates the causal role of the IL-33/ST2 pathway in dementia, with a focus on Alzheimer’s disease (AD), using Mendelian Randomisation (MR) methods. It seeks to understand whether the association between the IL-33/ST2 pathway and dementia is causal by leveraging genetic data from the UK Biobank.

### Background
Dementia, including Alzheimer's disease and vascular dementia, affects millions of people worldwide, with over 850,000 cases in the UK alone. Previous studies have implicated the immune-regulating IL-33/ST2 axis in the progression of Alzheimer's disease. This project aims to explore the genetic contribution of this axis to dementia, providing insights that may inform preventive or therapeutic measures.

### Methods
The project uses both individual-level and summary-level Mendelian Randomisation methods:
1. Individual-level MR: Utilises polygenic risk scores (PRS) as genetic proxies for IL-33 and ST2, combined with logistic regression to assess the impact on dementia outcomes.
2. Summary-level MR: Conducted using external genome-wide association studies (GWAS) to increase statistical power.
Key assumptions of Mendelian Randomisation, including relevance, exchangeability, and exclusion restriction, are considered to ensure the robustness of the analysis.

### Implementation
All relevant R code for data processing, MR analysis, and statistical modelling can be found in the [Scripts](Scripts) folder.

### Results
- ST2 PRS showed a significant association with its serum levels and with dementia outcomes in observational analysis, but Mendelian randomisation failed to establish a direct genetic link between ST2 or IL-33 PRS and dementia risk.
- Future studies may explore stratification by APOE-ε4 status or utilise tissue-specific protein data to overcome limitations in this study.

### Conclusion
The study highlights the relevance of the IL-33/ST2 axis in dementia, particularly in Alzheimer's disease. However, the results suggest no direct genetic determination of dementia risk via this pathway. Further research, including larger and more diverse datasets, is required to fully elucidate this relationship.

### Requirements
R
MendelianRandomization R package
UK Biobank data (restricted access)

### Acknowledgement
This project is a team collaboration as part of the Translational Data Science module (Teammates: Amin Moghaddam, Hanh Lan Bui, Lucie Frechin, Riya Nagar). It was carried out under the supervision of Dr. Verena Zuber and Dr. Roopal Desai.

Special thanks to Rin Wada and Joel Heller for their invaluable support in the technical aspects of the project.
