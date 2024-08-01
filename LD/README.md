# LD calculations between SNPs and SVs
### You will first have to concatenate VCF files that contain both the SNPs and SVs used in the GWAS analysis
```
bcftools concat \
--allow-overlaps \
--threads $SLURM_CPUS_PER_TASK \
-O z \
-o hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat.vcf.gz \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_sorted_gwas_SNPs.vcf.gz \
282_Peiffer2014Genetics_bimbam_w_header_20240716_bimbam_header.csv_bimbam_gwas_svs_reheader_sorted_hapmap_order.vcf.gz
```
```
bcftools sort \
-o hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat_sorted.vcf.gz \
-O z \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat.vcf.gz
```
```
tabix \
-p vcf \
-f \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat_sorted.vcf.gz
```
### Calculate LD between SVs and SNPs, making sure to first find SNPs that are found within SVs
```
bedtools intersect \
-a hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_sorted_gwas_SNPs.vcf.gz.bed \
-b 282_Peiffer2014Genetics_wo_header_bimbam_20240716.csv.bed \
-wa | cut -f4 > hapmap_gwas_SNPs_in_all_SVs.txt
```

```
plink \
--vcf hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_cat_sorted.vcf.gz \
--allow-extra-chr \
--exclude hapmap_gwas_SNPs_in_all_SVs.txt \
--ld-snp-list gwas_SVs.txt \
--ld-window 1000000 \
--ld-window-kb 1000 \
--ld-window-r2 0 \
--make-founders \
--out hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_ld \
--r2 gz dprime with-freqs \
--vcf-half-call m

```
### Convert LD file to make it easier to csv format for ease of handling
```
zcat hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_ld.ld.gz | sed -e "s~\s\+~,~g" | sed -e "s~^,~~g" | sed -e "s~,$~~g" > hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_ld.ld_sedcomma
```
### Split LD file into separate chromosomes
```
for f in {1..10} ; do head -1 hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_ld.ld_sedcomma > hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_chr${f}_ld.ld_sedcomma ; tail -n +2 hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_ld.ld_sedcomma | grep -P "^${f}," >> hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_chr${f}_ld.ld_sedcomma & done
```
### Run R script to find the closeset SNP for each SV (modified from https://github.com/HirschLabUMN/TE_variation/tree/master/TE-SNP_LD_methods)
```
for f in {1..10} ; do Rscript ld_R2_nsc_20240111.R hapmap321_gwas_svs_maf0.1_miss0.1_1mb_SNPs_in_all_SVs_gwas_SNPs_svs.vcf_chr${f}_ld.ld_sedcomma ld_distribution_results ${f} ; done
```

### Combine output to make master table of highest r$$^2$$ for each SV
```
head -1 plink_results_SNPs-highest-LD-TE_R2_chr3.ld > plink_results_SNPs-highest-LD-TE_R2_allchrs.ld ; for f in {1..10}; do tail -n +2 plink_results_SNPs-highest-LD-TE_R2_chr${f}.ld >> plink_results_SNPs-highest-LD-TE_R2_allchrs.ld; done
```
