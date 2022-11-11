# m.tuberculosis-research-code
Research code from article 'Feature selection and aggregation for antibiotic resistance GWAS in Mycobacterium tuberculosis: a comparative study'

## PFAM domains likelihood prediction



## Phenotype Heatmap

To get phenotype heatmap figure (Fig1) run jupyter notebook ./Pheno heatmap (Figure 1).ipynb

## Data Preprocessing

### Matrix deduplication

Before usage matrix gotten in the previous step it's needed to deduplicate feature by follow scripts in the data_preprocessing directory:

	python fix_duplicates.py {drug}
	python fix_test_duplicates.py {drug}

Specify input_dir and output_dir variables inside code script within configuration part to run these scripts

## Test Aggregated Features

Codes for the part are in the ./aggregated_features folder.

First of all, it's needed using logistic regression to get metrics for Indels & SNPs datasets by:

    python lr_snp_indel.py {drug}
    
For Indels & SNPs & aggregated by gene features datasets by:

    python lr_snp_indel_aggr.py {drug}
    
For Indels & SNPs & broken gene features datasets by:

    python lr_snp_broken.py {drug}
    
For Indels & SNPs & aggregated by PFAM features datasets by:

    python lr_snp_indel_pfam.py {drug}
    
For Indels & SNPs & without exculding mutations datasets by:

    python lr_snp_indel_thr1.py {drug}
    
And after executing all these scripts, run jupyter notebook 'Feature sets analysis.ipynb' to get TableS2-S10 (tables with comparation feature sets and dispersions for the sets)

## Feature Selection

### ABESS

Use this scripts in abess folder to run ABESS for each drug and each fold:

	Rscript run_abess.R {drug} {fold}

(Specify input_dir and output_dir variables inside code script within configuration part)

After that apply logistic regression on non-zeros selected features within abess folder:

	python abess_lr.py

(Specify dat_dir, result_dir and outptu_dir variables inside code script within configuration part to run these scripts)

And finally run jupyter notebook -- ABESS analysis -- to make Table5, TableS11, TableS17 (and some preparation for Table3). Change data_dir, data_dir_pre, output_abess_lr, final_output_abess if it's needed

### HHS

Before running HHS, it's needed to convert matrix to HHS format and make distance hamming files:

    python hamming.py {drug}
    python convert_to_hhs.py {drug}

After that, make run_multi_hhs.sh executive and run three time with different minimum frequency of mutations that lead to the resistance phenotype (f = 1, f = 3 and f = 5). Change it within configuration part within the file:

    ./run_multi_hhs.sh
    
Apply logstic regression on selected HHS features with different f (1, 3, 5) by:

    python hhs_lr.py
    
Finally, run jupyter notebook -- HHS analysis.ipynb -- to make TableS12, TableS18 (and some preparation for Table3)

All code is the ./hhs folder.

### Logistic regression with L1, SCAD regularization (pycasso package)

Codes for the part are in the ./lr folder.

First of all, run the below code with l1, scad and regulization ({method} = 'scad', 'mcp', 'l1') and gamma parameter (for scad and mcp {gamma} = 3, 5, 10, 15, and for l1 leave it empty):

    python py_picasso.py {drug} {method} {gamma}
    
And to make TableS13-15, TableS19-21 (and some preparation for Table3) run 'Pycasso analysis.ipynb'


### Logistic regression with elastic net regularization

Codes for the part are in the ./elastic_net folder.

Execute the below code with parameter a that equal 0.25, 0.5 and 0.75:

    python elastic_net.py {drug}
    
After that run jupyter notebook - Elastic net analysis.ipynb - to create TableS16, TableS22 (and some preparation for Table3)

### Feature metrics 

To generate Table 3 (The comparison of abilities of the feature selection methods to catch up features for the genes associated with known mechanisms of drug resistance) run jupyter notebook -- './Feature metrics (Table 3).ipynb'

### Supplement feature metrics

To generate Table S24 (The comparison of stability of feature selection by different method) run jupyter notebook -- './Supplement feature metrics (Table S24).ipynb'

## Population structure

### Making matrix with TreeBreaker features

Inside population_structure folder there are script to make matrix with TreeBreaker features to next usage ABESS feature selection:

	python make_db_with_str_features.py {drug}

After getting matrix, deduplicate them using instruction Data Preprocessing / Matrix deduplication and use script for selecting features by ABESS using instruction - Features Selection/ABESS. Specify input and output folders for all scripts!


