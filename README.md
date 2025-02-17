# binGO-GS

### Genomic Prediction of Arabidopsis thaliana Using SNP Subset Selected by Integrating GO-terms and Combinational Optimization in Bins
 
## Functional characteristics
We propose an improved GS method called binGO-GS, which integrates gene ontology (GO)-based prior biological knowledge with a novel bin-based variable selection algorithm to identify a subset of SNP markers that affect phenotypic variation, aiming to reduce genotyping costs and improve genomic prediction accuracy.Quantitative traits from two A. thaliana datasets with nearly 1000 samples and over half a million SNPs. Compared with using all the markers or randomly selected markers for prediction with support vector regression, the marker subset selected by binGO-GS achieved statistically significant improvements in prediction accuracy across all the traits, with p values of 0.0134 and 4.54×10-27, respectively.  For the six other reference models, binGO-GS also showed significant improvements in all six models compared with using all the markers in these models.  The markers selected by binGO-GS align with the polygenic hypothesis of minor-effect genes underlying complex quantitative traits, revealing that the selected markers for identical or similar morphological traits exhibit similar trends in quantity and distribution. 

## Installation procedure
### Getting Started
#### Prerequisites
```http
plink1.7
Rstudio
gtf file
R packages:
clusterProfiler
org.At.tair.db
data.table
rrBLUP
```
### Steps to Install
1. Clone the repository: 
```http
git clone https://github.com/ZhijunBioinf/binGO-GS.git
cd binGO-GS
Rscript Explore_num_SNPS.R
```
2. Install dependencies:
```http
pip install -r requirements.txt
```


## Usage
Organize your genotypic and phenotypic data and gtf annotation files:
```http
genotypic_train.csv
phenotypic_test.csv
```
### Train the Model
##### To train the model, run:
```http
Rscript Explore_num_SNPS.R
Rscript Find_target_subset.R
```
## Contributing

This project was developed by:


- Qingfang Ba (1468222359@qq.com) - Implementing
- Zhijun Dai (daizhijun@hunau.edu.cn) - Supervisor  
We welcome contributions from the community! Feel free to submit pull requests or raise issues.


