
## Alignments

### Create file containing SV-present and SV-absent alleles for each polymorphic SV within your ascertainment set
```
python anchorwave_include_sv_insertion_fastagen.py \
-t ../data/B73_Oh43_SV_summary_example.tsv \
-a B73 \
-b Oh43 \
-g Zm-B73-REFERENCE-NAM-5.0_chrs.fa \
-j Zm-Oh43-REFERENCE-NAM-1.0_chr.fa
```
#### Note: The above example uses a sample datatype for proof of concept.



### Download SRAs

### Trim fastq files using Trimmomatic
```
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -basein SRA_1.fastq -baseout SRA.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
```
##### Note: TruSeq3-PE-2.fa can be found [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa)

### Align with HISAT2

```

hisat2-build B73_Oh43_all_chrs_pseudoref.fasta B73_Oh43_all_chrs_pseudoref.fasta

hisat2 -q --phred33 --time --no-mixed --new-summary --no-discordant --summary-file ${SAMPLE_NAME}_summary_file.txt --met-stderr -k 1 -X $INSERT_SIZE -p $CPUS -x $HISAT2_INDEX --un ${SAMPLE_NAME}_un_unpaired_failed.fastq --un-conc ${SAMPLE_NAME}_unconc_failed_to_align_concordantly_%.fastq -1 ${FASTQ1} -2 ${FASTQ2} -S ${SAMPLE_NAME}.sam
```

## ReadID Sort Coordinate sort hisat2 output sam
```
samtools sort \
-@ $SLURM_CPUS_PER_TASK \
-n \
-o ${SAMPLE_NAME}_idsorted.sam \
${SAMPLE_NAME}.sam
```

## Fixmate to add info
```
samtools fixmate \
-@ $SLURM_CPUS_PER_TASK \
-r \
-m \
${SAMPLE_NAME}_idsorted.sam \
${SAMPLE_NAME}_idsorted_fixmate.sam
```

## Coordinate sort hisat2 output sam
```
samtools sort \
-@ $SLURM_CPUS_PER_TASK \
-o ${SAMPLE_NAME}_coordsorted_fixmate.sam \
${SAMPLE_NAME}_idsorted_fixmate.sam
```


## Remove PCR duplicates
```
samtools markdup \
-@ $SLURM_CPUS_PER_TASK \
-r \
-s \
-f ${SAMPLE}_markdup_stats.txt \
-u \
${SAMPLE_NAME}_coordsorted_fixmate.sam \
${SAMPLE_NAME}_coordsorted_fixmate_markdup.sam
```

## Filter to only include MAPQ score of 60
```
samtools view \
-q 60 \
-h \
-@ $SLURM_CPUS_PER_TASK \
-o ${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.sam \
${SAMPLE_NAME}_coordsorted_fixmate_markdup.sam
```
# SAM to BAM
```
samtools view \
-h \
-b \
-@ $CPUS \
-o ${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam \
${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.sam
```

# Index BAM
```
samtools index \
-b \
-@$SLURM_CPUS_PER_TASK \
${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam
```
# Check to make sure you are getting expected EOF
```
samtools quickcheck \
-v \
${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam
```

# Convert BAM to BED
```
bedtools bamtobed \
-split \
-i ${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam \
> ${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam_split.bed
```
# # bedtools coverage with left and right junction coordinates to get the number of reads that span over both the left and right junctions of SVs
```
bedtools coverage \
-a ../../anchorwave/anchorwave_pseudoref_generation/B73_Oh43_left_right_sv_insertionpoint_junctions.bed \
-b ${SAMPLE_NAME}_coordsorted_fixmate_markdup_mapq60.bam_split.bed \
> 33-16_bedtools_coverage_B73_Oh43_left_right_sv_insertionpoint_junctions.bed
```

## Create table for read depth at SV junctions for SV-present alleles and polymorphic site for SV-absent alleles 
```
python ../../../scripts/python_scripts/difference_reads_spanning_sv_junctions_insertion_point_NAM_genotypes_to_B73_Oh43_insertions.py \
../B73/B73_reads_align_B73_Oh43_all_alleles_pseudoref_k1_nodiscordant_normal_insert_size_coordsorted_fixmate_markdup_mapq60.bam_split_w_insertion_point.bed_bedtools_coverage_B73_Oh43_left_right_sv_insertionpoint_junctions.bed \
../Oh43/Oh43_reads_align_B73_Oh43_all_alleles_pseudoref_k1_nodiscordant_normal_insert_size_coordsorted_fixmate_markdup_mapq60.bam_split_w_insertion_point.bed_bedtools_coverage_B73_Oh43_left_right_sv_insertionpoint_junctions.bed \
33-16_bedtools_coverage_B73_Oh43_left_right_sv_insertionpoint_junctions.bed \
33-16 \
> 33-16_bedtools_coverage_B73_Oh43_left_right_sv_insertionpoint_junctions_comparison_table.bed
```
