# python3

'''
This is used for RIMSseq project.
Use to combine the mpileup counting generated by CountErrorMpileup.py in the same region listed in the annotation.

Require: bedtools executable in PATH

Logic:
(1) combine the Read1 and Read2 of positions included in the same region
Count the total and error for all the positions listed in input for the same region listed in the annotation file.
(2) with --context
top strand    ....TGCAcGGT...
bottom strand ....ACGTGcCA...
The lower case indicates CpG Cs on both strand. 
The CpG C on top strand is at pos=chr2:240523941
The CpG C on bottom strand is at pos=chr2:240523942

Read1 for 240523941 if Read1 mapped to the top strand: --REF C --BASE T --left 1 --right 0 return context AC
Read2 for 240523941 if Read1 mapped to the top strand: --REF G --BASE A --left 0 --right 1 return GT

Usage Example:
$python CountErrorRegion.py \
--Read1 Tar_NA12878_RIMseq_Rep1.R1.candidateStableRegion_C2T.CpG.count \
--Read2 Tar_NA12878_RIMseq_Rep1.R2.candidateStableRegion_G2A.CpG.count \
    --output Tar_NA12878_RIMseq_Rep1.candidateStableRegion.CpG.RegionCount \
--target PMID_25493099_Methylated_1kb.candidateStableRegion.gtf --context

Options:
--Read1: 
Read1 mpileup.count file generated using CountErrorMpileup.py, 1-coordination gft like format
--Read2: 
Read2 mpileup.count file generated using CountErrorMpileup.py, 1-coordination gft like format
Here I input Read1 and Read2 seperately because the context analysis for them is different.

--output: 1-coordination gft like format
The col1-col9 is the same as the annotation file.
The source and feature are the same as entries in annotation.
e.g. 
without --context col10 is the counting
chr1    ucsc    CpGisland       923591  926252  .       .       1-coordination  cpgNum=325;gcNum=1919;  count_total=46656;count_error=104;
with --context, col10 is the counting and col11 is the context
chr1    ucsc    CpGisland       923591  926252  .       .       1-coordination  cpgNum=325;gcNum=1919;  count_total=46656;count_error=104;  AC

--target: a gtf file containing the target regions for counting
e.g. CGI gtf annotation, 
/mnt/home/yan/Bo/RIMSeq/reference/ucsc_CpGisland.Blacklisted.mainchr.gtf
The counting/position intersecting with the same region in this annotation will be combined.
It the counting is outside the regions in this annotation, it will be ignored from counting.
This file should not contain a header line.

--context: Optional
Add this option to combine the counting in the context dependent manner, so the same region will have combined counting for different context.
To combine the Read1 and Read2 context based counting, the context of Read2 input is reverse complemented to match the Read1 context.
See CountErrorMpileup.py for details.
e.g.
Read1 context AC is cominbed with Read2 context GT for the same region.
'''
import os
import re
import argparse
import time
import datetime
from subprocess import check_call

__version__ = '2022.8.8'

def RC_seq(seq):
    dic={'A':'T', 'T':'A', 'G':'C','C':'G', 'N':'N'}
    seq = seq.upper()
    ls = list(seq)
    ls_new = []
    for letter in ls[::-1]:
        ls_new.append(dic[letter])
    return ''.join(ls_new)

class CountAnnotation:
    def __init__(self, Read1, Read2, output_file, annotation, context):
        self.output_file = output_file
        self.annotation = annotation
        self.Read1 = Read1
        self.Read2 = Read2
        self.context = context

        localtime = time.asctime(time.localtime())
        self.prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

        self.Run()

    def Run(self):

        # generate a tmp dir to save all the temp files, delete when finished
        cwd = os.getcwd()
        try:
            tmp = os.path.join(cwd, 'tmp{}'.format(self.prefix))
            print('Generate a tmp dir {}.'.format(tmp))
            os.mkdir(tmp)
            os.chdir(tmp)

            self.input_mpileupcount = self.CombineInput(self.Read1, self.Read2) # combined and sorted input file for intersection
            tempfile = self.Intersect()

            if not self.context: # Do not consider context for counting
                print("Peform context independent region counting.")
                self.Combine(tempfile)
            else: # context dependent counting
                print("Peform context dependent region counting." )
                print("The Read2 context will be reverse complemented first to combine with the Read1 context.")
                self.CombineContext(tempfile)
        except Exception as e: # print the error message
            print(e)
            print("Error! Error Counting failed.")
        finally: # clean up, delete the tmp dir when finished or failed
            os.chdir(cwd)
            check_call(['rm', '-r', '{}'.format(tmp)], shell=False)

    def extractTotal(self, word):
        return int(re.findall('count_total=(\d.*?);', word)[0])
    def extractError(self, word):
        return int(re.findall('count_error=(\d.*?);', word)[0])
    
    def CombineInput(self, Read1, Read2):
        '''
        with --context:
            combine all the Read1 files, combine all the Read2 files mean while reverse complement the context of Read2
        without --context:
            combine all the Read1 and Read2 files
        '''
        # sort the combined file for bedtools intersect
        if self.context:
            with open('combineInput.{}'.format(self.prefix), 'w') as output:
                with open(Read1) as f:
                    for line in f:
                        print(line.strip(), file=output, end='\n')
                # Reverse complement the context of Read2 for following analysis
                with open(Read2) as f:
                    for line in f:
                        line = line.strip().split('\t')
                        line[-1] = RC_seq(line[-1])
                        print('\t'.join(line), file=output, end='\n')

            command = 'sort -f -k1,1 -k4,4n -V {} > {}'.format('combineInput.{}'.format(self.prefix), 'combineInput.sorted.{}'.format(self.prefix))
            check_call(command, shell=True)
        else:
            with open('combineInput.{}'.format(self.prefix), 'w') as output:
                for input in [Read1, Read2]:
                    with open(input) as f:
                        for line in f:
                            print(line.strip(), file=output, end='\n')
            command = 'sort -f -k1,1 -k4,4n -V {} > {}'.format('combineInput.{}'.format(self.prefix), 'combineInput.sorted.{}'.format(self.prefix))
            check_call(command, shell=True) # with '>' can not use shell=False
          
        return 'combineInput.sorted.{}'.format(self.prefix)


    def Intersect(self):
        '''
        intersect region and mpileup file
        '''
        with open('bedtools_{}.sh'.format(self.prefix), 'w') as output:
            print('#!/bin/sh', file=output, end='\n')
            command = 'bedtools intersect -a {} -b {} -wa -wb > intersect.{}.temp'.format(self.annotation, self.input_mpileupcount, self.prefix)
            print(command, file=output, end='\n')
        
        try:
            check_call(['sh', 'bedtools_{}.sh'.format(self.prefix)], shell=False)
        except:
            print('Error when using bedtools intersect.')
        return 'intersect.{}.temp'.format(self.prefix)
        
    def Combine(self, intersectFile):
        '''
        combine the counts in the same region without considering context
        The region having no coverage count_total=0 will also be included in the output file.
        '''
        with open(self.annotation) as f:
            key = ['{}:{}-{}'.format(line.split()[0], line.split()[3], line.split()[4]) for line in f] # {(chr:start-end): count_total or count_error}
            dic_total = dict.fromkeys(key, 0)
            dic_error = dict.fromkeys(key, 0)
        with open(intersectFile) as f:
            for line in f:
                key = '{}:{}-{}'.format(line.split()[0], line.split()[3], line.split()[4])
                dic_total[key] += self.extractTotal(line.strip())
                dic_error[key] += self.extractError(line.strip())
        output = open(self.output_file, 'w')
        with open(self.annotation) as f:
            for line in f:
                line = line.strip().split('\t')
                key = '{}:{}-{}'.format(line[0], line[3], line[4])
                line.append('count_total={};count_error={};'.format(dic_total[key], dic_error[key]))
                print('\t'.join(line), file=output, end='\n')
        output.close()
        os.remove(intersectFile)
        return self.output_file

    def CombineContext(self, intersectFile):
        '''
        combine the counts in the same region considering context
        The region-context having no coverage will also be included in the output file.

        somehow bedtools intersect combine the counting and context, e.g.
        chr1    PMID_25493099   Ultra-stable    103553499       103554499       .       .       1-coordination  cg23707905;Methylated;  \
            chr1    Tar_NA12878_RIMseq_Rep1_MAPQ20.R1.Standard0     mpileup 103553886       103553886       .       -       \
                1-coordination  count_total=1;count_error=0;TC
        '''
        with open(intersectFile) as f:
            key = ['{}:{}-{};{}'.format(line.split()[0], line.split()[3], line.split()[4], line.strip().split()[-1].split(';')[-1]) for line in f] # chr1:103553499-103554499;TC
            dic_total, dic_error = dict.fromkeys(key, 0), dict.fromkeys(key, 0) # {(chr:start-end;context): count_total or count_error}
            context_combination = list(set([item.split(';')[1] for item in key]))
            f.seek(0)
            for line in f:
                key = '{}:{}-{};{}'.format(line.split()[0], line.split()[3], line.split()[4], line.strip().split()[-1].split(';')[-1])
                dic_total[key] += self.extractTotal(line.strip())
                dic_error[key] += self.extractError(line.strip())
        output = open(self.output_file, 'w')
        with open(self.annotation) as f:
            for line in f:
                line = line.strip().split('\t')
                ls_key = ['{}:{}-{};{}'.format(line[0], line[3], line[4], item) for item in context_combination]
                for key in ls_key:
                    temp = []
                    temp.extend(line)
                    temp.append('count_total={};count_error={};'.format(dic_total.get(key, 0), dic_error.get(key, 0)))
                    temp.append(key.split(';')[-1]) # context
                    print('\t'.join(temp), file=output, end='\n')
        output.close()
        os.remove(intersectFile)
        return self.output_file

def GetFilePath(input_file):
    '''
    return the abspath of the input_file
    if use relative path and input_file is not in cwd, the abspath will be wrong, return False
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
    parser.add_argument('--Read1', help='CountErrorMpileup.py output files for Read1', dest='Read1')
    parser.add_argument('--Read2', help='CountErrorMpileup.py output files for Read2', dest='Read2')
    parser.add_argument('--output', help='gtf file containing the counting results', dest='output_file')
    parser.add_argument('--target', help='gtf annotation file containing the target regions for error counting, e.g. CpG islands', dest='target')
    parser.add_argument('--context', help = 'add this option to perform context dependent counting', dest = 'context', action= 'store_true')
    
    args = parser.parse_args()
    return vars(args) # convert Namespace object into a dictionary

def main(dic_args):
    '''
    Here I use dic_args to faciliate running this function as an imported module instead of from command line
    dic_args {'Read1': , 'Read2': , 'output_file': , 'target': , 'context': True or False}

    convert to the full path of all the used files
    '''
    dic_args['Read1'] = GetFilePath(dic_args['Read1'])
    dic_args['Read2'] = GetFilePath(dic_args['Read2'])
    dic_args['output_file'] = os.path.abspath(dic_args['output_file'])
    dic_args['target'] = GetFilePath(dic_args['target'])

    print(datetime.datetime.now())
    print("The verion of CountErrorRegion.py script is: {}".format(__version__))
    print("The parameters used are:")
    print("Read1 count file: {}".format(dic_args['Read1']))
    print("Read2 count file: {}".format(dic_args['Read2']))
    print('Region annoation file: {}'.format(dic_args['target']))
    print("Output file: {}".format(dic_args['output_file']))
    
    CountAnnotation(dic_args['Read1'], dic_args['Read2'], dic_args['output_file'], dic_args['target'], dic_args['context']) # args.Read1: a list
    print('Done.')

    return 1

##-------mainbody
if __name__ == '__main__':
    main(Arg())

    
    
