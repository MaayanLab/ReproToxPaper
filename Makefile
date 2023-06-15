# run this to get python dependencies
.PHONY: pip-install
pip-install:
	pip install -r requirements.txt

# pure downloads

# more information about this: https://doi.org/10.1093/database/baab017
data/DrugRepurposingHub_moa_drugsetlibrary_name.dmt:
	curl -L https://maayanlab-public.s3.amazonaws.com/drugmonizome-dmts/DrugRepurposingHub_moa_drugsetlibrary_name.dmt -o $@

# more information about this: https://doi.org/10.1093/bioinformatics/bty060
data/Drugs_metadata.csv:
	curl -L https://maayanlab.cloud/L1000FWD/download/Drugs_metadata.csv -o $@

# more information about this: https://doi.org/10.1093/nar/gkac328
data/LINCS_small_molecules.tsv:
	curl -L https://s3.amazonaws.com/lincs-dcic/sigcom-lincs-metadata/LINCS_small_molecules.tsv -o $@

# more information about this: https://doi.org/10.1186%2Fs12859-022-04590-5
data/L1000_2021_drug_similarity.npz:
	curl -L https://appyters.maayanlab.cloud/storage/DrugShot/L1000_2021_drug_similarity.npz -o $@

# pre-processing

# TODO: these come from processing drugbank/chembl attrs
data/drugbank_drug_synonyms.dmt:
data/chembl_drugs.json:
data/2022-08-29-drug_synonmys.dmt:

# TODO: this comes from some rdkit processing
data/drug_attribute_mats.h5:

# TODO: this comes from sigcom lincs umap code
data/parametric-umap-all-meta.tsv:
data/lincs-umap-2021-08-24-chem-umap.tsv:

data/drugs-com.tsv:
	python 2022-09-09-drugs-com-crawl.py $@

# actual code for figures

.PHONY: 2023-04-25-id-mapping
2023-04-25-id-mapping: 2023-04-25-id-mapping.py input/pregnancy_category_D_and_X_v2.xlsx data/L1000_2021_drug_similarity.npz data/Drugs_metadata.csv data/LINCS_small_molecules.tsv data/drugbank_drug_synonyms.dmt data/chembl_drugs.json data/drugs-com.tsv
	python $<

data/2023-04-25-placenta-mapped.tsv: 2023-04-25-id-mapping
data/2023-04-25-D-X-mapped.tsv: 2023-04-25-id-mapping
data/2023-04-25-drugs-com.tsv: 2023-04-25-id-mapping

.PHONY: 2022-08-30-pregnancy-drug-preds
2022-08-30-pregnancy-drug-preds: 2022-08-30-pregnancy-drug-preds.py input/pregnancy_category_D_and_X_v2.xlsx input/Placenta_Barrier.xlsx data/LINCS_small_molecules.tsv data/parametric-umap-all-meta.tsv data/lincs-umap-2021-08-24-chem-umap.tsv data/DrugRepurposingHub_moa_drugsetlibrary_name.dmt
	python $<

.PHONY: 2023-04-28-benchmark
2023-04-28-benchmark: 2023-04-28-benchmark.py data/L1000_2021_drug_similarity.npz
	python $<

.PHONY: fig-2
fig-2: 2023-04-25-id-mapping

.PHONY: fig-3
fig-3: 2022-08-30-pregnancy-drug-preds

.PHONY: fig-4
fig-4: 2023-04-28-benchmark

# TODO
.PHONY: fig-s1
fig-s1:
