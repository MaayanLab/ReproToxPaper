# run this to get python dependencies
.PHONY: pip-install
pip-install:
	pip install -r requirements.txt


# provided

# this was extracted from ChEMBL 30 https://chembl.gitbook.io/chembl-interface-documentation/downloads
input/chembl_drugs.json:

# this file maps L1000 chemical signatures to UMAP-coordinates and is from https://maayanlab.cloud/sigcom-lincs/#/UMAPs
input/2021-08-30-chem-umap.tsv.gz:


# pure downloads

# more information about this: https://doi.org/10.1093/database/baab017
data/DrugRepurposingHub_moa_drugsetlibrary_name.dmt:
	curl -L https://maayanlab-public.s3.amazonaws.com/drugmonizome-dmts/DrugRepurposingHub_moa_drugsetlibrary_name.dmt -o $@

# more information about this: https://doi.org/10.1093/bioinformatics/bty060
data/Drugs_metadata.csv:
	curl -L https://maayanlab.cloud/L1000FWD/download/Drugs_metadata.csv -o $@

data/Probes_L1000_metadata.csv:
	curl -L https://maayanlab.cloud/l1000fwd/download/Probes_L1000_metadata.csv -o $@

# more information about this: https://doi.org/10.1093/nar/gkac328
data/LINCS_small_molecules.tsv:
	curl -L https://s3.amazonaws.com/lincs-dcic/sigcom-lincs-metadata/LINCS_small_molecules.tsv -o $@

data/cp_coeff_mat.gctx:
	curl -L https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gctx/cd-coefficient/cp_coeff_mat.gctx -o $@

# more information about this: https://doi.org/10.1186%2Fs12859-022-04590-5
data/L1000_2021_drug_similarity.npz:
	curl -L https://appyters.maayanlab.cloud/storage/DrugShot/L1000_2021_drug_similarity.npz -o $@

data/drugbank-full-database.xml:
	echo "You must download this from https://go.drugbank.com/releases/latest"

data/resources.dmt:
	curl -L https://s3.amazonaws.com/maayan-kg/reprotox/dmt/resources.dmt -o $@

# actual code for figures

.PHONY: 2021-08-30-ldp3-figures
2021-08-30-ldp3-figures: 2021-08-30-ldp3-figures.ipynb data/Probes_L1000_metadata.csv data/cp_coeff_mat.gctx
	python $<

# NOTE: this file is linked from `input` since it requires much computational resources
#  to produce but you can swap the two statements below to produce it instead
# data/2021-08-30-chem-umap.tsv.gz: 2021-08-30-ldp3-figures
data/2021-08-30-chem-umap.tsv.gz: input/2021-08-30-chem-umap.tsv.gz
	ln -s data/2021-08-30-chem-umap.tsv.gz input/2021-08-30-chem-umap.tsv.gz

.PHONY: 2022-09-09-drugs-com-crawl
2022-09-09-drugs-com-crawl: 2022-09-09-drugs-com-crawl.py
	python $<

2022-09-09-drugs-com-crawl: data/drugs-com.tsv

.PHONY: 2022-05-16-drugbank
2022-05-16-drugbank: 2022-05-16-drugbank.py data/drugbank-full-database.xml
	python $<

data/drugbank_enzyme.dmt: 2022-05-16-drugbank
data/drugbank_transporter.dmt: 2022-05-16-drugbank
data/drugbank_drug_synonyms.dmt: 2022-05-16-drugbank

.PHONY: 2022-05-16-dmt-prep-drugbank
2022-05-16-dmt-prep-drugbank: 2022-05-16-dmt-prep-drugbank.py data/drugbank_enzyme.dmt data/drugbank_transporter.dmt data/drugbank_drug_synonyms.dmt
	python $<

data/drugshot_autorif_enzyme_drugbank.tsv: 2022-05-16-dmt-prep-drugbank
data/drugshot_autorif_transporter_drugbank.tsv: 2022-05-16-dmt-prep-drugbank
data/drugshot_drugrif_enzyme_drugbank.tsv: 2022-05-16-dmt-prep-drugbank
data/drugshot_drugrif_transporter_drugbank.tsv: 2022-05-16-dmt-prep-drugbank

.PHONY: 2022-11-11-all-rdkit-features
2022-11-11-all-rdkit-features: 2022-11-11-all-rdkit-features.py data/L1000_2021_drug_similarity.npz data/LINCS_small_molecules.tsv data/drugbank_drug_synonyms.dmt data/drugbank_enzyme.dmt data/drugbank_transporter.dmt data/drugshot_autorif_enzyme_drugbank.tsv data/drugshot_autorif_transporter_drugbank.tsv data/drugshot_drugrif_enzyme_drugbank.tsv data/drugshot_drugrif_transporter_drugbank.tsv input/chembl_drugs.json
	python $<

.PHONY: 2022-11-28-more-features
2022-11-28-more-features: 2022-11-28-more-features.py 2022-11-11-all-rdkit-features
	python $<

data/drug_attribute_mats.h5: 2022-11-28-more-features
data/2022-08-29-drug_synonmys.dmt: 2022-11-28-more-features

.PHONY: 2023-04-25-id-mapping
2023-04-25-id-mapping: 2023-04-25-id-mapping.py input/pregnancy_category_D_and_X_v2.xlsx data/L1000_2021_drug_similarity.npz data/Drugs_metadata.csv data/LINCS_small_molecules.tsv data/drugbank_drug_synonyms.dmt input/chembl_drugs.json data/drugs-com.tsv
	python $<

data/2023-04-25-placenta-mapped.tsv: 2023-04-25-id-mapping
data/2023-04-25-D-X-mapped.tsv: 2023-04-25-id-mapping
data/2023-04-25-drugs-com.tsv: 2023-04-25-id-mapping

.PHONY: 2022-08-30-pregnancy-drug-preds
2022-08-30-pregnancy-drug-preds: 2022-08-30-pregnancy-drug-preds.py input/pregnancy_category_D_and_X_v2.xlsx input/Placenta_Barrier.xlsx data/LINCS_small_molecules.tsv data/2021-08-30-chem-umap.tsv.gz data/DrugRepurposingHub_moa_drugsetlibrary_name.dmt
	python $<

.PHONY: 2023-04-28-benchmark
2023-04-28-benchmark: 2023-04-28-benchmark.py data/L1000_2021_drug_similarity.npz
	python $<

2023-02-07-supervenn: 2023-02-07-supervenn.py data/resources.dmt
	python $<

.PHONY: fig-2
fig-2: 2023-04-25-id-mapping

.PHONY: fig-3
fig-3: 2022-08-30-pregnancy-drug-preds

.PHONY: fig-4
fig-4: 2023-04-28-benchmark

.PHONY: fig-s1
fig-s1: 2023-02-07-supervenn
