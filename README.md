# Metabolic modeling of Streptococcus pneumoniae


## Prepare S. Pneumonia strain T4 GEM models


### **Strategy 1: Modify existing GEM models** 

S. pneumoniae R6 model was obtained from [Dias (2019)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.01283/full)
```bash
!curl -L -o Dias2019_iDS372.xml "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000001090?filename=Dias2019_iDS372.xml"
```

Model: `data/Dias2019_iDS372.xml`


S. pneumoniae D39 model was obtained from [Pedram (2020)](https://pubmed.ncbi.nlm.nih.gov/32880642/) from Supplemntal Data `data/Additional_file_1.xlsx` and processed to `sbml` format using `scripts/prep_D39_model.py`

Model:  `data/iSPD39_model.xml`

TODO: 
- Map orthologs between strains
- Gap fill

### **Strategy 2: Build model from scratch from genome sequence**

T4 genome and annotation obtained from [ncbi](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/885/GCF_000006885.1_ASM688v1/)

TODO:
- Run CarveMe
