#!/bin/bash
##Job settings for without -t settings
#SBATCH --job-name=plink_LD_hapmap321_sv_vcf_svs_to_snps_1mb_miss0.1_gwas_SNPs_svs
#SBATCH -o plink_LD_hapmap321_sv_vcf_svs_to_snps_1mb_miss0.1_gwas_SNPs_svs.out
#SBATCH -e plink_LD_hapmap321_sv_vcf_svs_to_snps_1mb_miss0.1_gwas_SNPs_svs.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=catlinna@msu.edu

##Job Resources
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300gb
#SBATCH --time=7-00:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --reservation=centos_compute

cd $SLURM_SUBMIT_DIR

module purge

CONDA_BASE=$(conda info --base)

source $CONDA_BASE/etc/profile.d/conda.sh

conda activate plink


#############
# VARIABLES #
#############


##############
## COMMANDS ##
##############


# this is now looking at all sv ids
# ld-window is kind of perfect becaus that is just over the higehst snp count for a  chrom
plink \
--vcf hapmap321_all_chrs_peiffer_intersection_genos_picard_liftovervcf_v3_to_v4_to_v5NAM_w_gwas_SNPs_gwas_svs_cat_sorted_20240722.vcf.gz \
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

# Below is from https://github.com/HirschLabUMN/TE_variation/blob/master/TE-SNP_LD_methods/scripts/plink_te-snp.sh
#plink \
#--file analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.plk \
#--make-founders \
#--r2 gz dprime with-freqs \
#--ld-window-r2 0 \
#--ld-window 1000000 \
#--ld-window-kb 1000 \
#--out analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime \
#--allow-extra-chr
