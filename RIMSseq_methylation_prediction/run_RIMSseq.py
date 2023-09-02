# -*- coding: utf-8 -*-
# python3

'''
Use to create a bash script to predict the methylation for certain regions based on RIMSseq results.
Here perform context dependent calibaration for the methylation prediction for human genomic regions, e.g. CpG islands.

Requirement:
samtools, bedtools exacutable in $PATH
need regex, pandas, numpy and statistics modules

Usage:
Run run_RIMSseq.py to generate a bash script containing all the steps for human regional methylation prediction using RIMSseq.
Suggest running the generated bash script on a cluster since it may take a while to finish the counting.

Steps contained:
(1) remove Read1-Read2 overlapping, optional
(2) split the Read1 and Read2 and perform mpileup counting without --no-BAQ.
(3) count error for each position using CountErrorMpileup.py: Use --REF C --BASE T for Read1, --REF G --BASE A for Read2.
(4) combine the error counting based on regions using CountErrorRegion.py in the context-dependent manner: the 0-standard, 100-standard and target region (e.g. CGI). 
(5) Based on the error rate in 0-standard and 100-standard, predict the methylation level in target region using PredictMethylation.py in the context-dependent manner.

If starting with bam file with removal of PCR duplicates, performing steps (1)-(5)
$python run_RIMSseq.py --input name.dedup_reads.bam --removeOverlap \
        --name name --MAPQ 10 --vcf mnt/home/yan/Bo/AccurateSeq/NA12878_WGS/variant_bowtie2/HG001_WGS_gatktmp/HG001_WGS_filtered.vcf \
		--RegionAnnotation ~/annotation/ucsc_CpGisland.Blacklisted.mainchr.gtf \
        		--CpGAnnotation ~/annotation/ucsc_CpGisland.CpGsite.mainchr.gtf \
        			--WGBS ~/annotation/ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt

If starting with bam file with removal of PCR duplicates and removal of Read1-Read2 overlapping, performing steps (2)-(5)
$python run_RIMSseq.py --input /mnt/home/yan/Bo/RIMSeq/220318_RIMSeq/bwamem/name.deover_reads.bam \
    --name name --MAPQ 0 \
        --vcf /mnt/home/yan/Bo/AccurateSeq/NA12878_WGS/variant_bowtie2/HG001_WGS_gatktmp/HG001_WGS_filtered.vcf \
            --RegionAnnotation /mnt/home/yan/Bo/RIMSeq/reference/ucsc_CpGisland.Blacklisted.mainchr.gtf \
                --CpGAnnotation /mnt/home/yan/Bo/RIMSeq/reference/ucsc_CpGisland.CpGsite.mainchr.gtf

If Starting with Read1 and Read2 mpileup files, performing steps (3)-(5)
$python run_RIMSseq.py --R1mpileup name.R1.mpileup --R2mpileup name.R1.mpileup \
    --name name \
    --vcf /mnt/home/yan/Bo/AccurateSeq/NA12878_WGS/variant_bowtie2/HG001_WGS_gatktmp/HG001_WGS_filtered.vcf \
        --RegionAnnotation /mnt/home/yan/Bo/RIMSeq/reference/ucsc_CpGisland.Blacklisted.mainchr.gtf \
            --CpGAnnotation /mnt/home/yan/Bo/RIMSeq/reference/ucsc_CpGisland.CpGsite.mainchr.gtf

Options:

--input:
a bam file containing both Read1 and Read2; can be:
(1) bam file having PCR duplicates removed and Read1-Read2 overlapping removed; this bam file is used as input for mpileup of Read1 and Read2;
(2) bam file having PCR duplicates removed; need to add --removeOverlap; sort by name will be performed in the removal of overlapping step.

--R1mpileup, --R2mpileup: 
If Read1 and Read2 mpileup files are both provided, skip the mpileup() step.
Use the provided mpileup files to count the error and predict methylation.

--removeOverlap: Optional
Add this option to perform the Read1 Read2 overlapping removal step,
generate a bam file name.deover_reads.bam in the current working directory.

--name: name of the generated output files
e.g. name
The corresponding output files generated are:
name.R1.mpileup and name.R2.mpileup
name.CpG_MethylationPrediction

--MAPQ: INT, default 10
MPAQ cutoff used in mpileup --min-MQ.

--vcf: Default None
a vcf file used to remove variant positions from error counting
With default value, this is the command for CountErrorMpileup.py: --vcf None

--RegionAnnotation:
The gtf file containing target genomic regions for methylation prediction, used for CountErrorRegion.py --target option.
Only The CpG sites in the annotated region are used for error counting.
e.g. human CGI or promoter region annoation file

--CpGAnnotation: 
The gtf file containing CpG sites in the target regions defined by --RegionAnnotation option.
The CpG sites of this file is used for CountErrorMpileup.py --CpG option.
This file could be generated using ExtractCpGSite.py with --RegionAnnotation file as input file.

--cutoff: INT, Default 3000
Total coverage cutoff for RIMSseq used for PredictMethylation_Context.py --cutoff option.
Only the region having total>=cutoff will be saved in the output file CpG_MethylationPrediction.

--R1cycle or --R2cycle:
cycles will be ignored from error counting by CountErrorMpileup.py --cycle option.

Note:
(1) The input ultra stable region annotation files for standard0 and standard100 provided in ~/annotation folder, 
    which are defined in createScript().
(2) The temporary files that are generated:
name.R1.bam, name.R2.bam
name.R1.Standard100.C2T_sitecount, name.R2.Standard100.G2A_sitecount, name.Standard100.regioncount
name.R1.Standard0.C2T_sitecount, name.R2.Standard0.G2A_sitecount, name.Standard0.regioncount
name.R1.AnnotatedRegion.C2T_sitecount, name.R2.AnnotatedRegion.G2A_sitecount, name.AnnotatedRegion.regioncount
(3) Do not run multiple jobs using the same name defined by --name option, to avoid overwriting the temporary files.
The temporary files generated during samtools sort are named by both --name and the running time to avoid conflicts.
'''

import argparse
import os
import time

script_dir = os.path.join(os.path.dirname( __file__ ), 'src')
annotation_dir = os.path.join(os.path.dirname( __file__ ), 'annotation')

header = '''#!/bin/sh'''

__version__ = '2022.06.19' # version of the script

class CountError:
    '''
    Use to count the error (C2T for R1 and G2A for R2) based on RIMSseq results, mpileup counting for R1 and R2
    The command is appended to the end of script
    bamfile: Removal of overlapping between Read1 and Read2

    Output having error count of give regions, '{}.{}.regioncount'.format(self.name, self.source), 
    e.g. name.Standard100.regioncount, name.Standard0.regioncount, name.CGI.regioncount
    '''
    def __init__(self, bamfile, reference, name, source, script, MAPQ, BQ, vcf, CpG, region, percent_cutoff, count_cutoff, R1cycle, R2cycle):
        self.reference = reference
        self.name = name # e.g. name_CGI
        self.source = source # e.g. Standard100, Standard0, CGI
        self.bamfile = bamfile
        self.script = script

        self.MAPQ, self.BQ = MAPQ, BQ
        self.region = region # e.g. PMID_25493099_Methylated_1kb.candidateStableRegion.gtf
        self.vcf, self.CpG = vcf, CpG # CpG e.g. PMID_25493099_Methylated_1kb.candidateStableRegion.CpGsite.gtf
        self.percent_cutoff, self.count_cutoff = percent_cutoff, count_cutoff # parameters for CountErrorMpileup.py

        self.R1cycle = R1cycle # a list containing cycles ignored from error counting
        self.R2cycle = R2cycle # a list containing cycles ignored from error counting

    def mpileup(self, R1mpileup, R2mpileup):
        '''
        Select R1 and R2 for mplieup
        Output files:
        name.R1.mpileup
        name.R2.mpileup
        '''
        if R1mpileup and R2mpileup:
            return R1mpileup, R2mpileup
        else: # mppileup files are not provided, generated using samtools
            R1, R2 = '{}.R1.bam'.format(self.name), '{}.R2.bam'.format(self.name)
            with open(self.script, 'a') as output:
                command = 'samtools view -h -f 64 {} -b > {}'.format(self.bamfile, R1)
                print(command, file=output, end='\n')
                command = 'samtools view -h -F 64 {} -b > {}'.format(self.bamfile, R2)
                print(command, file=output, end='\n')
                # mpileup file name, e.g. name_MAPQ20.R1.mpileup
                command = 'samtools mpileup -f {} {} --min-MQ {} --min-BQ {} --output-BP-5 --output {} \
--no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends'.format(self.reference, R1, self.MAPQ, self.BQ, '{}_MAPQ{}.R1.mpileup'.format(self.name, self.MAPQ))
                print(command, file=output, end='\n')
                command = 'samtools mpileup -f {} {} --min-MQ {} --min-BQ {} --output-BP-5 --output {} \
--no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends'.format(self.reference, R2, self.MAPQ, self.BQ, '{}_MAPQ{}.R2.mpileup'.format(self.name, self.MAPQ))
                print(command, file=output, end='\n')
    
                command = 'rm {}'.format(R1)
                print(command, file=output, end='\n')
                command = 'rm {}'.format(R2)
                print(command, file=output, end='\n')

                print(file=output, end='\n')
            return '{}_MAPQ{}.R1.mpileup'.format(self.name, self.MAPQ), '{}_MAPQ{}.R2.mpileup'.format(self.name, self.MAPQ)

    def CountMpileup(self, R1mpileup, R2mpileup):
        '''
        counting region with --context
        output file e.g. name_CGI.R1.Standard0.C2T_sitecount
        '''
        with open(self.script, 'a') as output:
            # R1 CtoT
            if self.R1cycle:
                command = 'python {}/CountErrorMpileup.py --input {} --output {} \
--reference {} --left 1 --right 0 --REF C --BASE T \
--name {} --vcf {} --CpG {} --percentage_cutoff {} --minimum_count {} --cycle {}'.format(script_dir, R1mpileup, '{}.R1.{}.C2T_sitecount'.format(self.name, self.source), \
    self.reference, '{}.R1.{}'.format(self.name, self.source), self.vcf, self.CpG, self.percent_cutoff, self.count_cutoff, ' '.join(self.R1cycle))
                print(command, file=output, end='\n')
            else:
                command = 'python {}/CountErrorMpileup.py --input {} --output {} \
--reference {} --left 1 --right 0 --REF C --BASE T \
--name {} --vcf {} --CpG {} --percentage_cutoff {} --minimum_count {}'.format(script_dir, R1mpileup, '{}.R1.{}.C2T_sitecount'.format(self.name, self.source), \
    self.reference, '{}.R1.{}'.format(self.name, self.source), self.vcf, self.CpG, self.percent_cutoff, self.count_cutoff)
                print(command, file=output, end='\n')
            # R2 GtoA
            if self.R2cycle:
                command = 'python {}/CountErrorMpileup.py --input {} --output {} \
--reference {} --left 0 --right 1 --REF G --BASE A \
--name {} --vcf {} --CpG {} --percentage_cutoff {} --minimum_count {} --cycle {}'.format(script_dir, R2mpileup, '{}.R2.{}.G2A_sitecount'.format(self.name, self.source), \
    self.reference, '{}.R2.{}'.format(self.name, self.source), self.vcf, self.CpG, self.percent_cutoff, self.count_cutoff, ' '.join(self.R2cycle))
                print(command, file=output, end='\n')
            else:
                command = 'python {}/CountErrorMpileup.py --input {} --output {} \
--reference {} --left 0 --right 1 --REF G --BASE A \
--name {} --vcf {} --CpG {} --percentage_cutoff {} --minimum_count {}'.format(script_dir, R2mpileup, '{}.R2.{}.G2A_sitecount'.format(self.name, self.source), \
    self.reference, '{}.R2.{}'.format(self.name, self.source), self.vcf, self.CpG, self.percent_cutoff, self.count_cutoff)
                print(command, file=output, end='\n')
            
            command = 'python {}/CountErrorRegion.py --Read1 {} --Read2 {} --output {} \
--target {} --context'.format(script_dir, '{}.R1.{}.C2T_sitecount'.format(self.name, self.source), \
    '{}.R2.{}.G2A_sitecount'.format(self.name, self.source), \
        '{}.{}.regioncount'.format(self.name, self.source), self.region)
            print(command, file=output, end='\n')

            # delete tempory files used for error counting
            for item in ['{}.R1.{}.C2T_sitecount'.format(self.name, self.source), '{}.R2.{}.G2A_sitecount'.format(self.name, self.source)]:
                print('rm {}'.format(item), file=output, end='\n')

            print(file=output, end='\n')

        return '{}.{}.regioncount'.format(self.name, self.source)


def PredictMethylation(name, Standard100_regioncount, Standard0_regioncount, CGI_regioncount, WGBS, cutoff, script):
    '''
    Predict methylation based on error counting of 100% standard, 0% standard and CGI region

    cutoff: coverage cutoff for methylation prediction of given regions
    '''
    with open(script, 'a') as output:
        command = 'python {}/PredictMethylation_Context.py --target {} --standard-0 {} --standard-100 {} \
--output {}.CpG_MethylationPrediction --WGBSannotation {} --cutoff {} \
--stdev'.format(script_dir, CGI_regioncount, Standard0_regioncount, Standard100_regioncount, name, WGBS, cutoff)
        print(command, file=output, end='\n')

        # delete tempory files used for error counting
        for item in [Standard100_regioncount, Standard0_regioncount, CGI_regioncount]:
            print('rm {}'.format(item), file=output, end='\n')

    return '{}.CpG_MethylationPrediction'.format(name)

def createScript(args):
    '''
    create a script to run the pipeline for RIMseq methylation prediction
    '''
    reference = args.reference
    vcf = args.vcf

    localtime = time.asctime(time.localtime())
    prefix = '{}_{}'.format(name, ''.join(localtime.split()[-2].split(':'))) # e.g. Tar_NA12878_RIMseq_Rep2_CGI_151542

    # Region annotation: 100%-methylation-standard-Region, 0%-methylation-standard-Region, UCSC-CGI or other regions
    region_list = ['{}/PMID_25493099_Methylated_1kb.candidateStableRegion.gtf'.format(annotation_dir), \
                   '{}/PMID_25493099_Methylated_1kb.candidateStableRegion.gtf'.format(annotation_dir), \
                    args.regionFile]

    # CpG site annotation: 100%-methylation-standard-CpGsites, 0%-methylation-standard-CpGsites, UCSC-CGI-CpGsites or CpG sites of other regions
    CpG_list = ['{}/PMID_25493099_Methylated_1kb.candidateStableRegion.CpGsite.gtf'.format(annotation_dir), \
                '{}/PMID_25493099_Methylated_1kb.candidateStableRegion.nonCpGsite.gtf'.format(annotation_dir), \
                    args.CpGFile]

    source_list = ['Standard100', 'Standard0', 'AnnotatedRegion']

    MAPQ, BQ = args.MAPQ, 30
    percent_cutoff, count_cutoff = 0.5, 10
        
    name = args.name # 'Tar_NA12878_RIMseq_Rep2_CGI'
    script = '{}.MethylPrediction.script'.format(name)

    with open(script, 'w') as output:
        print(header, file=output, end='\n')

        if args.overlap:
            print('echo Remove the Read1 and Read2 overlapping.', file=output, end='\n')
            tempsam = os.path.basename(args.bamfile).replace('.bam', '.sortbyname.sam') # e.g. name.dedup_reads.sortbyname.sam
            command = 'samtools sort -n -T {} -o {} {}'.format(prefix, tempsam, args.bamfile) # the temporary files generated during sorting: prefix.nnnn.bam
            print(command, file=output, end='\n')
            command = 'python {}/TrimOverlappingReadPair.py --input {} --output {}'.format(script_dir, tempsam, tempsam.replace('.sam', '.deover_reads.sortbyname.sam'))
            print(command, file=output, end='\n')
            
            bamfile = ('{}.deover_reads.bam'.format(name)) # name.deover_reads.bam overlapping removed, sorted bam file used for error counting
            command = 'samtools sort -o {} {}'.format(bamfile, tempsam.replace('.sam', '.deover_reads.sortbyname.sam'))
            print(command, file=output, end='\n')
            command = 'rm {}'.format(tempsam)
            print(command, file=output, end='\n')
            command = 'rm {}'.format(tempsam.replace('.sam', '.deover_reads.sortbyname.sam'))
            print(command, file=output, end='\n')
        else:
            bamfile = args.bamfile # '/mnt/home/yan/Bo/RIMSeq/220318_RIMSeq/bowtie2/mpileup/Tar_NA12878_RIMseq_Rep2.dedup_deover_reads.bam'
    
    with open(script, 'a') as output:
            # check whether can find the input bam file or not
            if not args.R1mpileup and not args.R2mpileup:
                if not bamfile: # check bam file is provided or not
                    print('Need to provide a bamfile as input.')
                    quit()
                else: # check bam file exists or not
                    command = 'if [ -f \"{}\" ]; then true; else echo input not found && exit 1; fi'.format(bamfile)
                    print(command, file=output, end='\n')
                    print(file=output, end='\n')
            
            print('date', file=output, end='\n')
            print('echo =========', file=output, end='\n')
            if not args.R1mpileup and not args.R2mpileup:
                print('echo Context dependent Prediction methylation of Annotated Regions based on bam file:', file=output, end='\n')
                print('echo {}'.format(bamfile), file=output, end='\n')
                print('echo MAPQ, BQ used: {}, {}'.format(MAPQ, BQ), file=output, end='\n')
            else:
                print('echo Context dependent Prediction methylation of Annotated Regions based on mpileup files:', file=output, end='\n')
                print('echo {} and {}'.format(args.R1mpileup, args.R2mpileup), file=output, end='\n')

            print('echo vcf file used:', file=output, end='\n')
            print('echo {}'.format(vcf), file=output, end='\n')
            print('echo Region Annotation file used:', file=output, end='\n')
            print('echo {}'.format(args.regionFile), file=output, end='\n')
            print('echo CpG sites Annotation file used:', file=output, end='\n')
            print('echo {}'.format(args.CpGFile), file=output, end='\n')
            print('echo RIMSseq methylation prediction file generated is: {}'.format('{}.CpG_MethylationPrediction'.format(name)), file=output, end='\n')
            if args.R1cycle:
                print('echo The following illumina sequencing cycles will be ignored from Read1 error counting: {}'.format(' '.join(args.R1cycle)), file=output, end='\n')
            if args.R2cycle:
                print('echo The following illumina sequencing cycles will be ignored from Read1 error counting: {}'.format(' '.join(args.R2cycle)), file=output, end='\n')
            print('echo =========', file=output, end='\n')
            print(file=output, end='\n')

    output_list = [] # Standard100_regioncount, Standard0_regioncount, CGI_regioncount for PredictMethylation function
    Flag = True
    for region, CpG, source in zip(region_list, CpG_list, source_list):
        temp = CountError(bamfile, reference, name, source, script, MAPQ, BQ, vcf, CpG, region, percent_cutoff, count_cutoff, args.R1cycle, args.R2cycle)
        if Flag: # only generate mpileup file once
            R1mpileup, R2mpileup = temp.mpileup(args.R1mpileup, args.R2mpileup)
            Flag = False
        output_list.append(temp.CountMpileup(R1mpileup, R2mpileup))

        with open(script, 'a') as output:
            print(file=output, end='\n')
    PredictMethylation(name, output_list[0], output_list[1], output_list[2], args.WGBS, args.cutoff, script)

    with open(script, 'a') as output:
        print(file=output, end='\n')
        print('date', file=output, end='\n')

    return script

##-------mainbody
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='bam file used for error counting', dest='bamfile', required=False)
    parser.add_argument('--R1mpileup', help='R1 mpileup file', dest='R1mpileup', default=None)
    parser.add_argument('--R2mpileup', help='R2 mpileup file', dest='R2mpileup', default=None)
    parser.add_argument('--removeOverlap', help = 'remove R1 and R2 overlapping', dest = 'overlap', action= 'store_true')


    parser.add_argument('--name', help='name for prediction file', dest='name')
    parser.add_argument('--MAPQ', help='MAPQ cutoff for mpileup', dest='MAPQ', type=int, default=10)
    parser.add_argument('--cutoff', help='total coverage cutoff for RIMSseq error counting', dest='cutoff', type=int, default=3000)
    
    parser.add_argument('--reference', help='reference genome fa file used for mapping', dest='reference')

    parser.add_argument('--vcf', help='vcf file used to remove SNP sites from counting', dest='vcf', default=None)
    parser.add_argument('--RegionAnnotation', help='gtf annotation of regions used for methylation prediction', dest='regionFile')
    parser.add_argument('--CpGAnnotation', help='gtf annotation of CpG sites contained in the annotation regions', dest='CpGFile')
    parser.add_argument('--WGBS', help='WGBS file used for standard_100 calibration, e.g. ENCODE_human_WGBS_ultraStable_1kb_methylation.summary.txt', dest='WGBS')

    parser.add_argument('--R1cycle', nargs="+", help='sequencing cycles that are removed from R1 error couting', dest='R1cycle', default=None)
    parser.add_argument('--R2cycle', nargs="+", help='sequencing cycles that are removed from R2 error couting', dest='R2cycle', default=None) 

    args = parser.parse_args()

    createScript(args)


