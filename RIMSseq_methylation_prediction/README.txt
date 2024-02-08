================
RIMS-seq analysis
Author: Bo Yan
Date: Aug 24, 2023
================

This tool is used to predict the methylation levels of human genomic regions, such as CpG islands, through the application of RIMS-seq.

# Principle of RIMS-seq and methylation prediction

RIMS-seq deaminates 5-Methylcytosince (5mC) to thymine (T) and cytosine (C) to uracil with different conversion rate.
The converted uracil will be removed using USER enzyme from Illumina sequencing, therefore the CtoT error reflects the methylation level.
We have shown that the CtoT error rate and methylation level have a linear relationship based on linear regression analysis.
The deamination rate of 5mC is about 1%, so for human genome sequencing, we do not expect an accurate methylation quantification at base resolution due to the sequencing depth limitation.
Instead we can accurately predict the regional methylation level (e.g. CpG islands methylation level) based on the regional CtoT error rate.

Because the deamination rate of 5mC varies between each library preparation and sequencing, we use ultra stable regions in the human genome as internal controls to calibrate the linear relationship.
The benefit of using internal control is that it has the same level of background error rate. 
In addition, it works for exome targeted enrich sequencing.  

# Analysis pipeline to use RIMS-seq to extract methylation and variant information for human samples, including

(1) Trim the Illumina adapter sequence
(2) Map the read pairs to human hg38 genome using bowtie2 default alignment with the addition of RG header for variant calling
(3) Remove the unmapped reads and PCR duplicates using Picard tools
(4) Call germline variant using gatk
(5) Predict regional methylation using RIMSseq_methylation_prediction

# Annotation files provided in ~/annotation folder

(1) ucsc_CpGisland.Blacklisted.mainchr.gtf: human hg38 CpG islands from UCSC browser on autosomes and chrX with blacklist regions removed, gtf format
(2) ucsc_CpGisland.CpGsite.mainchr.gtf: CpG sites in the CpG islands, gtf format
(3) PMID_25493099_Methylated_1kb.candidateStableRegion.gtf: ultra stable regions of hg38 human genome
(4) PMID_25493099_Methylated_1kb.candidateStableRegion.CpGsite.gtf: CpG sites in the ultra stable regions, used as fully methylated control
(5) PMID_25493099_Methylated_1kb.candidateStableRegion.nonCpGsite.gtf: non-CpG cytosine sites in the ultra stable regions, used as unmethylated control
(6) ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt: the methylation level of human ultra stable regions based on ENCODE WGBS data


# Python scripts provided

RIMSseq_methylation_prediction/run_RIMSseq.py

RIMSseq_methylation_prediction/src:
(1) TrimOverlappingReadPair.py
(2) CountErrorMpileup.py
(3) CountErrorRegion.py
(4) PredictMethylation_Context.py

# Requirement

All the python scripts are based on python 3;
Need python3 module regex, pandas, numpy and statistics;
Need samtools and bedtools executable in PATH.


# Usage Example

^ Steps contained in the methylation analysis:
(1) remove Read1-Read2 overlapping, optional
(2) split the Read1 and Read2 and perform mpileup counting without --no-BAQ.
(3) count error for each position using CountErrorMpileup.py: Use --REF C --BASE T for Read1, --REF G --BASE A for Read2.
(4) combine the error counting based on regions using CountErrorRegion.py in the context-dependent manner: the 0-standard, 100-standard and target region (e.g. CGI). 
(5) Based on the error rate in 0-standard and 100-standard, predict the methylation level in target region using PredictMethylation.py in the context-dependent manner.


^ Execute the following command for methylation prediction:

If starting with bam file with removal of PCR duplicates, performing steps (1)-(5)
$python run_RIMSseq.py --input name.dedup_reads.bam --removeOverlap \
    --name name --MAPQ 10 --reference hg38.fa \
            --RegionAnnotation ~/annotation/ucsc_CpGisland.Blacklisted.mainchr.gtf \
                --CpGAnnotation ~/annotation/ucsc_CpGisland.CpGsite.mainchr.gtf \
			--WGBS ~/annotation/ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt

If starting with bam file with removal of PCR duplicates and removal of Read1-Read2 overlapping, performing steps (2)-(5)
$python run_RIMSseq.py --input /mnt/home/yan/Bo/RIMSeq/220318_RIMSeq/bwamem/name.deover_reads.bam \
    --name name --MAPQ 0 --reference hg38.fa \
            --RegionAnnotation ~/annotation/ucsc_CpGisland.Blacklisted.mainchr.gtf \
                --CpGAnnotation ~/annotation/ucsc_CpGisland.CpGsite.mainchr.gtf \
			--WGBS ~/annotation/ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt

If Starting with Read1 and Read2 mpileup files, performing steps (3)-(5)
$python run_RIMSseq.py --R1mpileup name.R1.mpileup --R2mpileup name.R1.mpileup \
    --name name --reference hg38.fa \
            --RegionAnnotation ~/annotation/ucsc_CpGisland.Blacklisted.mainchr.gtf \
                --CpGAnnotation ~/annotation/ucsc_CpGisland.CpGsite.mainchr.gtf \
			--WGBS ~/annotation/ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt

Note:
Execute the provided command to generate a bash script named NAME.MethylPrediction.script, which will contain all the necessary commands for methylation prediction.

Subsequently, execute the generated bash file either locally or on the cluster. It is recommended to run this on the cluster, as the process may take a significant amount of time to complete.

Modify the generated script for:
	Adding proper header lines as required for running on the cluster.
	Select the appropriate Python environment with the required modules installed by activating the virtual environment (e.g., using 'conda activate') or by modifying the Python path.


^ Options:

--input: a bam file
Normally use a Bam file with removal of PCR duplicates

--name: name for the bash script and the output files
A bash script NAME.MethylPrediction.script is generated containing all the commands required for the methylation prediction.

These files are generated during the analysis:
NAME.deover_reads.bam, a bam file with removal of Read1 and Read2 overlapping
NAME_MAPQ10.R1.mpileup, Read1 samtools mpileup file
NAME_MAPQ10.R2.mpileup, Read2 samtools mpileup file

NAME.CpG_MethylationPrediction, human genomic region methylation prediction, e.g.
<region><total><error><Methylation><stdevMethylation>
chr1-10638243-10640041  29781   299     0.7952499339199948      0.04579085321566869
Region: 1-coordination
Methylation: 0-1 range
stdevMethylation: standard deviation of methylation prediction.

--MAPQ: Int, Default value 10
MAPQ cutoff for samtools mpileup counting

--vcf: Optional, a vcf file
If provided, variant positions in this vcf file are ignored from error counting.
Can use the variant calling results generated by gatk.

--RegionAnnotation: a gtf file
A gtf file containing the target genomic regions for methylation prediction, 
e.g. human CpG islands annotation is provided as ~/annotation/ucsc_CpGisland.Blacklisted.mainchr.gtf
This gtf file should contain 9 columns having:
Column1 chromosome
Column3 1-coordination start of region
Column4 1-coordination end of region
The output is the prediction of methylation level of these genomic regions.

--CpGAnnotation: a gtf file
A gtf file containing the CpG sites in the genomic regions listed in the --RegionAnnotation file.
Only the positions listed in this --CpGAnnotation file are used for error counting.
e.g. human CpG sites in CpG islands is provided as ~/annotation/ucsc_CpGisland.CpGsite.mainchr.gtf

--reference: a fasta file
human hg38 reference genome that is used for mapping.
The contigs in the --RegionAnnotation file should match this genome file.

--WGBS: Optional
File containing the methylation evaluation of the ultra stable regions by other benchmark methods.
The methylation based on ENCODE human WGBS is provided as ~/annotation/ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt
If provided, this file is used for calibrating the error rate of 100% methylation standard.

--cutoff: Int, Default 3000
Total coverage cutoff for RIMS-seq used as PredictMethylation_Context.py --cutoff.
Only the region having total>=cutoff will be saved in the output file CpG_MethylationPrediction.

--R1cycle or --R2cycle: Optional, using space as delimiter
If provided, these cycles will be ignored from error counting by CountErrorMpileup.py --cycle.
e.g. --cycle 0 1
The Illumina sequencing circle is 0-coordination, so 0 corresponds to the first cycle for both read1 and read2 in the fastq files that are used for mapping.
