# python3

'''
This is used for RIMSseq project.
Use to parse the mpileup file to count the REFtoBASE error (e.g. CtoT error on read) in the context (read context) dependent manner.

Require: bedtools executable in PATH

Usage Example:
$python CountErrorMpileup.py --input Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.mpileup \
    --output Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.Standard100.C2T_sitecount \
    --reference GCA_000001405.15_GRCh38_no_alt_spike.fa \
    --left 1 --right 0 --REF C --BASE T --name Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.Standard100 \
    --vcf ~/WGS_NA12878_RIMseq_Rep1_gatktmp/WGS_NA12878_RIMseq_Rep1_filtered.vcf \
    --CpG PMID_25493099_Methylated_1kb.candidateStableRegion.CpGsite.gtf \
    --percentage_cutoff 0.5 --minimum_count 10

Logic:
use --REF C --BASE T for RIMSseq deamination counting:
    error means the number of reads having T at a position that should be C in the read (fastq file); 
    total means the number of reads having T or C at a position that should be C, which is different from coverage on this position (also containing A or G).

Options:
--input: mpileup file
This mpileup file needs to be generated by samtools mpileup with the following option: 
    --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends --output-BP-5
Its format is shown as below:
    <chr><pos 1-coordination><col3 Reference=top genome><coverage><Read Bases><Base quality in ASCII characters><1-based positions within the read><context>

--name: Not required
The name used in source column in output.
If not provided, use the input file name as source.

--reference: reference genome used for mapping
fasta file for bedtools getfasta, the same file that is used for mapping

--left, --right: int bp, Default 0 for both to shut down the context calling
upstream or downstream of the postion, used for bedtools slop -l and -r which returns the seq of (chr, pos-l, pos+b)
l=1 and r=1 returning a 3bp context that the target base is in the middle

--REF, --BASE: type of error in the read, Default C, T
Default means that the base should be C based on the reference genome but it is T in the fastq.
Need to choose different setting for Read1 and Read2 for differnt purpose.
e.g. RIMseq deamination counting use Read 1 --REF C -- BASE T and Read2 --REF G --BASE A.

--output: output file
The output file has 1-coordination gft like format, containing the error counting.
(1) The positions having no coverage (count_total=0) are not shown in the output.
(2) source is decided by --name; feature is mpileup
(3) context is in the last column

If --left or --right !=0, call context
chr2    Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.Standard100   mpileup 240523629       240523629       .       +       1-coordination  count_total=14;count_error=0;   GC
The context in the last column is based on the strand to reflect the context in read for context dependent error counting.
--left 0 --right 0 will shut down the context calling
chr2    Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.Standard100   mpileup 240523629       240523629       .       +       1-coordination  count_total=14;count_error=0;


Options controlling the error counting:
--vcf: a vcf file for known SNP, Default None
If provided, varaint positions with column7 FILTER=PASS in the vcf file will not be used for error counting.
    Add --noVcfFilter to remove all the positions listed in the vcf file from error counting, in regardless of column 7 FILTER status.

--CpG: a gtf file containing the targeted C positions (e.g. CpG sites), Default None
If provided, only use/report the positions provided in this gtf file for error counting. 
For CpG sites, the annotation file should contain the positions on both strands for the CpG pair.

--cycle: 0-index sequencing cycles that will not be used for error counting, Default None
e.g. --cycle 73 74
If provided, remove bases at certain cycle(s) from error counting. 
e.g. --cycle 74, means remove the last sequencing cycle (cycle=74) for 75bp sequencing trimmed reads for both read1 and read2.
This corresponds to 75 in mpileup file <1-based positions within the read> column.

--percentage_cutoff: float, --minimum_count: Int
Default --percentage_cutoff 0.5 --minimum_count 10
If count_total>= minimum_count and error rate (count_error/count_total) at a position>=percentage_cutoff, counting/result at this position is removed from output.
Add this option to remove the weird loci having too many error, which could be due to SNP or mapping.
By default most non SNP positions will not be removed.
Add --minimum_count option to avoid removing the positions having low coverage but by chance having an error 
Set --percentage_cutoff 1.1 will shut down this filter in regardless of the --minimum_count value.

Note:
(1) running time estimation
For 217M bam file having 6.4 million of R1, it takes less than 5min to finish 'samtools mpileup', 
and it takes about 5min to count the mpile up file (about 2.1G).
(2) The counting in this script can be called both from command line and as an imported module as following:
    import CountErrorMpileup
    CountErrorMpileup.RunCountError(dic_args)
Here dic_args is a dictionary saving the require parameters, 
e.g. {'input_file': Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.mpileup, 'output_file': ,'reference': , 'l': 1, 'r': 1,  ...}.
(3) Prepare full path for all the files, or file without any path which means file is under cwd, do not use relative path e.g. './input_file'.
(4) Generate a tmp dir to save all the temp files, delete the tmp dir at the end.
'''
import os
import time
import re
from subprocess import check_call
import argparse

__version__ = '2023.06.19' # version of the script

class CountError:
    '''
    Functions used to count the REFtoBASE error and get the context.
    For Mpileup REF (column 3) = REF, context is based on top genome (5->3' direction controled by -l and -r in getfasta -s)
    For Mpileup REF (column 3) = RC of REF, context is based on the bottom genome (5->3' direction controled by -l and -r in getfasta -s)
    '''
    def __init__(self, input_file, REF, BASE, output_file, reference, r, l, cycle, exp, vcf, CpG, percentage_cutoff, minimum_count, noVcfFilter):
        # all the files have full path after conversion in RunCountError(0)
        self.input_file, self.output_file = input_file, output_file
        self.cycle = cycle
        self.REF, self.BASE = REF, BASE
        self.reference, self.r, self.l = reference, r, l
        self.vcf, self.CpG = vcf, CpG
        self.percentage_cutoff, self.minimum_count = percentage_cutoff, minimum_count
        self.exp = exp # source in gtf
        self.noVcfFilter = noVcfFilter # True or False

        localtime = time.asctime(time.localtime())
        self.prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

        self.Run()

    def Run(self):
        '''Run the class'''
        
        # Quit if the output file already exists
        if os.path.exists(self.output_file):
            print("The output file {} already exists. Quit".format(self.output_file))
            quit()
        
        # generate a tmp dir to save all the temp files, delete when finished
        cwd = os.getcwd()
        try:
            tmp = os.path.join(cwd, 'tmp{}'.format(self.prefix))
            print('Generate a tmp dir {}.'.format(tmp))
            os.mkdir(tmp)
            os.chdir(tmp)

            tempgtf = self.contionalSelectionForOther()
            if self.r>0 or self.l>0: # Add context in the last column
                print("Add context information in the last column.")
                print("Left and Right for context extraction: {} and {}".format(self.l, self.r))
                self.getFasta(tempgtf, self.output_file)
            else: # Do not have context added 
                print("Do not add context.")
                os.rename(tempgtf, self.output_file)
        except Exception as e: # print the error message
            print(e)
            print("Error! Error Counting failed.")
        finally: # clean up, delete the tmp dir when finished or failed
            os.chdir(cwd)
            check_call(['rm', '-r', '{}'.format(tmp)], shell=False)

        return self.output_file

    def filter(self, col_base, col_readpos, ls_cycle):
        '''
        col_base: ...,,,,, no separator
        col_readpos: '63,61,57,54,53,53,53...', separator is ','
        ls_cycle: a list containing the illumina sequencing cycles that should be removed from counting

        sequencing cycle filter, remove the read bases for certain cycles if ls_cycle is not empty
        cycle 0-index, 1-based positions 1-index. So 1-based positions = cycle + 1
        return a list containing the read bases with unwanted base at certain cycles removed

        soft clip and cycle cycle:
        See ~/Bo_tools_README/SAMFormat.txt section 1.5 for explanation.
        '''
        ls = []
        ls_readpos = col_readpos.split(',')
        for i in range(0, len(ls_readpos)):
            if str(int(ls_readpos[i])-1) not in ls_cycle:
                ls.append(col_base[i])
        return ''.join(ls)

    def count(self, ls_base, col3, REF, BASE):
        '''
        col3 is the column3 in mpileup, which is the reference of top genome of this position
        count the CtoT error following the logic
        '''
        dic_RC = {'C': 'G', 'G': 'C', 'A':'T', 'T':'A', 'N':'N'}
        if col3 == REF: # REF=C
            count_error = len(re.findall(BASE, ls_base)) # BASE=T
            count_total = count_error + len(re.findall('\.', ls_base)) # BASE=C
            return (count_total, count_error, '+')
        elif col3 == dic_RC[REF]: # dic_RC[REF]=G
            count_error = len(re.findall(dic_RC[BASE].lower(), ls_base)) # dic_RC[BASE].lower()=a
            count_total = count_error + len(re.findall(',', ls_base))
            return (count_total, count_error, '-')
        else:
            return (-1, -1, -1)

    def contionalSelectionForOther(self):
        '''
        Select counting from certain illumina sequencing cycles;
        Select positions in targeted regions (e.g. CGI) for counting;
        Remove positions on SNP site from counting.
        Report positions having counting above percentage_cutoff and minimum_count.
        '''
        if self.vcf:
            with open(self.vcf) as f:
                if self.noVcfFilter: # Use all the positions in VCF file
                    key = [(line.strip().split()[0], line.strip().split()[1]) for line in f if not line.startswith('#')]
                else: # Use positions with FILTER=PASS in VCF file
                    key = [(line.strip().split()[0], line.strip().split()[1]) for line in f if not line.startswith('#') and line.split('\t')[6]=='PASS']
                    print('There are {} positions in vcf file that are removed from error counting.'.format(len(key)))
                self.dic_variant = dict.fromkeys(key) # {(chr, 1-coordination position)}, python3.7 Dict keeps insertion order

        if self.CpG:
            with open(self.CpG) as f:
                # since it is 1bp position, do not consider strand
                key = [(line.strip().split()[0], line.strip().split()[3]) for line in f if not line.startswith('#')]
                self.dic_CpG = dict.fromkeys(key, (0, 0)) # {(chr, 1-coordination position): (count_total at this pos, count_error at this pos)}
                if len(self.dic_CpG.keys())==0: # in case the annotation file only has header line
                    print("There are no targeted sites in the annotation file.")
                    quit()

        output = open('{}.gtf'.format(self.prefix), 'w')
        with open(self.input_file) as f: 
            for line in f:
                line = line.strip().split('\t')
                
                if line[3] == '0': # coverage <1
                    continue
                if self.CpG:
                    if not (line[0], line[1]) in self.dic_CpG:
                        continue
                if self.vcf:
                    if (line[0], line[1]) in self.dic_variant:
                        continue
                
                if self.cycle:
                    filtered_base = self.filter(line[4], line[6], self.cycle)
                else:
                    filtered_base = line[4]

                count_total, count_error, strand = self.count(filtered_base, line[2], self.REF, self.BASE)
                if count_total>0: # Do not save positions having no coverage
                    if count_total<self.minimum_count or float(count_error)/count_total<self.percentage_cutoff:
                        temp = [line[0], self.exp, 'mpileup', line[1], line[1], '.', strand, '1-coordination', 'count_total={};count_error={};'.format(count_total, count_error)]
                        print('\t'.join(temp), file=output, end='\n')
        output.close()
        return '{}.gtf'.format(self.prefix)

    def getFasta(self, input_gtf, output_file):
        '''
        Add context based on --left or --right.
        If the requested number of bases exceeds the boundaries of the chromosome, bedtools slop will clip the feature accordingly.
        '''
        fai = self.reference + '.fai'
        if not os.path.exists(fai):
            print('Fai file is not found.')
            quit()

        with open('bedtools_{}.sh'.format(self.prefix), 'w') as output:
            print('#!/bin/sh', file=output, end='\n')
            command = 'bedtools slop -i {} -g {} -l {} -r {} -s | \
bedtools getfasta -fi {} -bed - -bedOut -s > {}'.format(input_gtf, fai, self.l, self.r, self.reference, '{}.fa'.format(self.prefix))
            print(command, file=output, end='\n')
        try:
            check_call(['sh', 'bedtools_{}.sh'.format(self.prefix)], shell=False)
        except:
            print('Error when using bedtools getfasta.')
        
        # paste the context in upper case into the gtf
        output = open(output_file, 'w')
        with open(input_gtf) as f1, open('{}.fa'.format(self.prefix)) as f2:
            line1 = f1.readline().strip()
            line2 = f2.readline().strip()
            while line1:
                temp = line1.split('\t')
                temp.append(line2.split('\t')[-1].upper())
                print('\t'.join(temp), file=output, end='\n')
                line1 = f1.readline().strip()
                line2 = f2.readline().strip()
        output.close()
            
        return output_file

def ConvertNone(x):
    '''
    x is string None/none -> return boolean None
    x is a boolean None -> return x as boolean None
    Does not work for args.cycle which is a list
    '''
    if x: # x is string or list
        if isinstance(x, str):
            x = None if x.lower() == 'none' else x
    return x

def GetFilePath(input_file):
    '''
    return the abspath of the input_file
    if use relative path and input_file is not in cwd, the abspath will be wrong, return False
    Do not use this function for output file since it is not created yet
    '''
    filepath = os.path.abspath(input_file)
    
    if os.path.exists(filepath):
        return filepath
    else:
        print('Can not find file: {}'.format(filepath))
        return None

def Arg():
    '''parser if run from command line'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='mpileup file', dest='input_file')
    parser.add_argument('--name', help='used as source in output file', dest='name', required=False)
    parser.add_argument('--output', help='output gtf file', dest='output_file')

    parser.add_argument('--reference', help='reference genome fa file used for mapping', dest='reference')

    parser.add_argument('--vcf', help='vcf file containing positions that are removed from error counting', dest='vcf', default=None)
    parser.add_argument('--CpG', help='gtf file containing positions that are reported for error counting', dest='CpG', default=None)
    parser.add_argument('--cycle', nargs="+", help='illumina sequencing cycles that are removed from error couting, 0-index', dest='cycle', default=None) 
    
    parser.add_argument('--percentage_cutoff', help='positions having coverage>=minimum_count and error_rate>=percentage_cutoff are discarded', \
                        dest='percentage_cutoff', default=0.5, type=float)
    parser.add_argument('--minimum_count', help='positions having coverage>=minimum_count and error_rate>=percentage_cutoff are discarded', \
                        dest='minimum_count', default=10, type=int)
    
    parser.add_argument('--REF', help='The theoretical base based on the reference genome', dest='REF', default='C')
    parser.add_argument('--BASE', help='The real based in the fastq file', dest='BASE', default='T')
    
    parser.add_argument('--left', help='upstream bp subtracted from the start coordinate for bedtools slop', dest='l', default=0, type=int)
    parser.add_argument('--right', help='downstream bp added to the end coordinate for bedtools slop', dest='r', default=0, type=int)

    parser.add_argument('--noVcfFilter', help='add this to remove all the positions in vcf instead of only FILTER=PASS positions for error counting', dest='noVcfFilter', \
        action='store_true', required=False)
    
    args = parser.parse_args()
    return vars(args) # convert Namespace object into a dictionary

def main(dic_args):
    '''
    Provide information and Run CountError.
    Here I use dic_args to faciliate running this function as an imported module instead of from command line
    dic_args {'input_file': Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.mpileup, 'output_file': ,'reference': , 'l': 1, 'r': 1, 'name': name}

    convert to the full path of all the used files
    '''
    # print(datetime.datetime.now())
    print("The verion of CountErrorMpileup.py script is: {}".format(__version__))
    print("The parameters used are:")

    dic_args['CpG'] = ConvertNone(dic_args['CpG'])
    dic_args['vcf'] = ConvertNone(dic_args['vcf'])

    dic_args['input_file'] = GetFilePath(dic_args['input_file'])
    dic_args['output_file'] = os.path.abspath(dic_args['output_file'])
    dic_args['reference'] = GetFilePath(dic_args['reference'])
    if not dic_args['input_file'] or not dic_args['reference']: # file does not exist or path is not correct
        quit()

    print("Input mpileup file: {}".format(dic_args['input_file']))
    print("Output file: {}".format(dic_args['output_file']))
    print("percentage_cutoff and minimum_count: {} and {}".format(dic_args['percentage_cutoff'], dic_args['minimum_count']))
    print("REF and BASE: {} and {}".format(dic_args['REF'], dic_args['BASE']))

    if 'noVcfFilter' in dic_args:
        noVcfFilter = True
    else:
        noVcfFilter = False
    if dic_args['vcf']:
        vcf = GetFilePath(dic_args['vcf'])
        print("vcf file containing positions that are ignored for error counting: {}".format(vcf))
        if not noVcfFilter:
            print('Only positions with FILTER=PASS in vcf file are ignored.')
    else:
        vcf = None
        print('No vcf file provided containing positions that are ignored for error counting.')
    
    if dic_args['CpG']:
        CpG = GetFilePath(dic_args['CpG'])
        print("File containing positions that are reported for error counting: {}".format(CpG))
    else:
        CpG = None
    
    # if name not provided, use the suffix of input file as name
    exp = dic_args.get('name', os.path.basename(dic_args['input_file']).split('.')[0])
    # check illumina cycle that will be ignored
    if dic_args['cycle']:
        print("The following illumina sequencing cycles will be ignored for error counting: {}".format(','.join(dic_args['cycle'])))
    
    print()

    CountError(dic_args['input_file'], dic_args['REF'], dic_args['BASE'], dic_args['output_file'], dic_args['reference'], \
dic_args['r'], dic_args['l'], dic_args['cycle'], exp, vcf, CpG, dic_args['percentage_cutoff'], dic_args['minimum_count'], noVcfFilter)

    print("Done.\n")

    return 1

##-------mainbody
if __name__ == '__main__':
    # run from command line
    main(Arg())


    

