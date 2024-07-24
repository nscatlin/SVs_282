## GEMMA associations

### Load GEMMA
```
conda env create -f gemma.yml
conda activate gemma
```

### Create kinship matrix
```
gemma -g 282_Peiffer2014Genetics_wo_header_bimbam.csv -p 282_Peiffer2014Genetics_phenotype_wo_extra_fields.tsv -o 282_Peiffer2014Genetics_phenotype_genotype_kinship_matrix_maf0.1_miss0.1 -gk 2 -miss 0.1 -maf 0.1 -outdir 0.1
```

### Run associations for each phenotype
```
for f in {1..11}; do gemma -g 282_Peiffer2014Genetics_wo_header_bimbam_20240716.csv -p 282_Peiffer2014Genetics_phenotype_wo_extra_fields_20240716.tsv -k ./0.1/282_Peiffer2014Genetics_phenotype_genotype_kinship_matrix_maf0.1_miss0.1.sXX.txt -a gemma_annotation_nodups_sorted.txt -lmm 4 -n ${f} -maf 0.1 -o ${f}_all_lmm_options_miss0.1 -miss 0.1 -outdir 0.1 > ${f}_gemma_assoc.out 2> ${f}_gemma_assoc.err & done
```

### Find significant associations
```
for f in {1..11}; do sed -e "s~PHENONUM~${f}~g" peiffer_man_qq_plot_miss0.1_pheno_tmp.R > peiffer_man_qq_plot_miss0.1_pheno_${f}.R ; Rscript peiffer_man_qq_plot_miss0.1_pheno_${f}.R > peiffer_man_qq_plot_miss0.1_pheno_${f}.out 2> peiffer_man_qq_plot_miss0.1_pheno_${f}.err & ; done
```
