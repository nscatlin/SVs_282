# LD calculations between SNPs and SVs
### You will first have to concatenate VCF files that contain both the SNPs and SVs used in the GWAS analysis
```
bcftools concat \
--allow-overlaps \
--threads $SLURM_CPUS_PER_TASK \
-O z \
-o hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat_20240716.vcf.gz \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_sorted_gwas_SNPs.vcf.gz \
282_Peiffer2014Genetics_bimbam_w_header_20240716_bimbam_header.csv_bimbam_gwas_svs_reheader_sorted_hapmap_order.vcf.gz
```
```
bcftools sort \
-o hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat_sorted_20240716.vcf.gz \
-O z \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat.vcf.gz
```
```
tabix \
-p vcf \
-f \
hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_concat_sorted_20240716.vcf.gz
```

```
