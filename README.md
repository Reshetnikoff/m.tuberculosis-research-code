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


## Population structure

### Making matrix with TreeBreaker features

Inside population_structure folder there are script to make matrix with TreeBreaker features to next usage ABESS feature selection:

	python make_db_with_str_features.py {drug}

Specify input_dir, output_dir, split_dir, feature_dir and pheno_dir variables inside code script within configuration part to run these scripts


