# m.tuberculosis-research-code
Research code from article 'Feature selection and aggregation for antibiotic resistance GWAS in Mycobacterium tuberculosis: a comparative study'

## PFAM domains likelihood prediction





## Data preprocessing

### Matrix deduplication

Before usage matrix gotten in the previous step it's needed to deduplicate feature by follow scripts in the data_preprocessing directory:

	python fix_duplicates.py {drug}
	python fix_test_duplicates.py {drug}

Specify input_dir and output_dir variables inside code script within configuration part to run these scripts

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

All code is inside ./hhs folder.

### Logistic regression with L1, SCAD regularization (pycasso package)

Codes for the part are in the ./lr folder.

First of all, run below code with l1, scad and regulization ({method} = 'scad', 'mcp', 'l1') and gamma parameter (for scad and mcp {gamma} = 3, 5, 10, 15, and for l1 leave it empty):

    python py_picasso.py {drug} {method} {gamma}
    
And to make TableS13-15, TableS19-21 (and some preparation for Table3) run 'Pycasso analysis.ipynb'

## Population structure

### Making matrix with TreeBreaker features

Inside population_structure folder there are script to make matrix with TreeBreaker features to next usage ABESS feature selection:

	python make_db_with_str_features.py {drug}

Specify input_dir, output_dir, split_dir, feature_dir and pheno_dir variables inside code script within configuration part to run these scripts


