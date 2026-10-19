#! /bin/bash -l
# # L segment
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_L.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --output LASV_GENBANK_20250725_curated.L.partitioned.beast

# # S segment
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_S.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --output LASV_GENBANK_20250725_curated.S.partitioned.beast

# # L segment no nigeria
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_L.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --skip-nigeria \
#     --output LASV_GENBANK_20250725_curated.L.noNIG.beast

# # S segment no nigeria
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_S.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --skip-nigeria \
#     --output LASV_GENBANK_20250725_curated.S.noNIG.beast

# L segment
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_L.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --segment L \
#     --output LASV_GENBANK_20250725_curated.L.partitioned.beast2

# # S segment
# python preparexml.py \
#     --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_S.renamed.msa \
#     --template xml-templates/template_codon_partitioned.xml \
#     --segment S \
#     --output LASV_GENBANK_20250725_curated.S.partitioned.beast2

# L segment no nigeria
conda run -n nf-core python preparexml.py \
    --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_L.renamed.msa \
    --template xml-templates/template_codon_partitioned.xml \
    --segment L \
    --skip-countries NIG NGA CIV BEN TGO MLI \
    --output runs/L/LASV_GENBANK_20250725_curated.L.noNIG-CIV.beast2

# S segment no nigeria
conda run -n nf-core python preparexml.py \
    --alignment ../../data/alignments/final/LASV_GENBANK_20250725_curated.added_background_S.renamed.msa \
    --template xml-templates/template_codon_partitioned.xml \
    --segment S \
    --skip-countries NIG NGA CIV BEN TGO MLI \
    --output runs/S/LASV_GENBANK_20250725_curated.S.noNIG-CIV.beast2