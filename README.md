# RIMS-seq2

**RIMS-seq2** can identify the human genomic variants and predict the regional methylation level (e.g. CpG islands methylation). <br>

Analysis pipeline to use RIMS-seq to extract methylation and variant information for human samples: <br>
(1) Trim the Illumina adapter sequence <br>
(2) Map the read pairs to human hg38 genome using bowtie2 default alignment with the addition of RG header for variant calling <br>
(3) Remove the unmapped reads and PCR duplicates using Picard tools <br>
(4) Call germline variant using gatk <br>
(5) Predict regional methylation using RIMSseq_methylation_prediction <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;See RIMSseq_methylation_prediction/README.txt for instructions for use. <br>
