# python 3

'''
Use to trim the overlapping region if Read1 and Read2 overlap.

logic:
If gtf_end of forward mapping read >= gtf_start of reverse mapping read, the pair of reads overlap.
    For overlapping read pair:
        Trim the 3'end of the forward mapping read to remove the overlapping region;
        Trim SEQ, QUAL, adjust the cigar and MD tag accrodingly;
        Trim the 3'end soft clip from the forward mapping if exist but Leave the 5'end soft clip in the forward mapping;
        The trimmed reap pairs have custom XO tag: XO:Z:matetrim or XO:Z:readtrim.
If no overlap, keep orignial Read1, Read2

Usage Example:
$python TrimOverlappingReadPair.py --input sortbyname.sam --output deover.sam --removeEmptyMapping

--input: a sam file
    the sam file containing paired mapping, sorted by name.
    cigar should only contain 'MDIS', otherwise may raise error.
--output: a sam file
    the sam file has trimmed read pairs.
    The trimmed reap pairs have custom XO tag: XO:Z:matetrim or XO:Z:readtrim.
    Output has the same order as input.
    Output has empyt entries (mapping having no matching) that may raise error for bedtools, see Note(1) for details.
--removeEmptyMapping: Optional
    Add this option to remove the mapping having no matching, which raises error for bedtools, see Note(1) for details.
    This is not necessary for samtools mpileup but necessary for bedtools.
    With this option, the reads may not paired in the output file.

Note:
Term used in this script:
InsertLength of a mapping: the length of the mapping part to the reference, equals to:
    1-coordination start-end+1;
    cigar sum(M|D);
    Also corresponds to the MDtag.
SEQLength: the length of fastq seq in Col10 SEQ, equals to:
    cigar sum(M|I|S)

(1) See BedTools_README.txt 1.8 section for examples and explaination.
If Read1 and Read2 have the same start in other words complete overlapping, generate one mapping having '*' in col10.
before trimming, read pair in /mnt/home/yan/Bo/CbDa/Novaseq_071823/bwamem/AcDa01_TET_3h_2ul_05.dedup_reads.bam:
A00336:A00336:H7H2GDRX3:1:1104:12328:12790	99	lambda	1	60	61M	=	1	56	GGGTGGTGACCTCGTGGGTTTTTGCTATTTATGAAAATTTTCTGGTTTAAGGTGTTTCTGT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	MD:Z:3C2C7C7C19C9C5C2	PG:Z:MarkDuplicates	RG:Z:A00336:H7H2GDRX3:1	NM:i:7	AS:i:31	XS:i:21
A00336:A00336:H7H2GDRX3:1:1104:12328:12790	147	lambda	1	60	56M	=	1	-56	GGGTGGTGACCTCGTGGGTTTTTGCTATTTATGAAAATTTTCTGGTTTAAGGTGTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	MD:Z:3C2C7C7C19C9C3	PG:Z:MarkDuplicates	RG:Z:A00336:H7H2GDRX3:1	NM:i:6	AS:i:30	XS:i:21
after trimming:
A00336:A00336:H7H2GDRX3:1:1104:12328:12790	99	lambda	1	60	0M	=	1	56	*	*	MD:Z:0	PG:Z:MarkDuplicates	RG:Z:A00336:H7H2GDRX3:1	NM:i:7	AS:i:31	XS:i:21	XO:readtrim
A00336:A00336:H7H2GDRX3:1:1104:12328:12790	147	lambda	1	60	56M	=	1	-56	GGGTGGTGACCTCGTGGGTTTTTGCTATTTATGAAAATTTTCTGGTTTAAGGTGTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	MD:Z:3C2C7C7C19C9C3	PG:Z:MarkDuplicates	RG:Z:A00336:H7H2GDRX3:1	NM:i:6	AS:i:30	XS:i:21	XO:matetrim
The trimmed forward mapping Flag 99 with '0M' will be removed from output when --removeEmptyMapping is present.

(2) Add the custom header, e.g.
    @PG	ID:TrimOverlappingReadPair	VN:20230730	CL:python TrimOverlappingReadPair.py --input AcDa01_TET_3h_2ul_05.dedup_reads.sortbyname.sam --output AcDa01_TET_3h_2ul_05.dedup_reads.deover.sam
(3) Add the custom XO tag in the trimmed read pairs:
    XO:Z:readtrim or XO:Z:matetrim, note need to use 'XO:Z:' format, otherwise raising error by samtools view.
'''
try:
    import regex
    import re
    import numpy as np
    from collections import deque
    import argparse
    import time
    import os
except:
    print('module error')
    quit()

__version__ = '2023.07.30'

def DetectSort(input_sam):
    '''
    Detect whether the sam file is sorted by name using header
    '''
    with open(input_sam) as f:
        if f.readline().strip().split()[-1] != 'SO:queryname':
            print("The input sam file is not sorted by name.")
            return False
        else:
            return True

class TrimOverlapping:
    def __init__(self, input_sam, output_sam, remove0S):

        localtime = time.asctime(time.localtime())
        self.prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

        self.input_sam = input_sam
        self.output_sam = output_sam

        # Run automatically when create class object
        temp_sam = self.GoThroughSam(self.input_sam, 'deover_{}.sam'.format(self.prefix))
        if remove0S:
            print("Remove the empty mapping from sam file.")
            self.RemoveEmptyMapping(temp_sam, self.output_sam)
            os.remove(temp_sam)
        else:
            os.rename(temp_sam, self.output_sam)

    def Flag16(self, Flag):
        '''
        return True if Flag contains Flag 16=2^4, Flase if not
        '''
        if len(str(np.binary_repr(int(Flag))))<5: # Flag < 16
            return False
        else:
            return str(np.binary_repr(int(Flag)))[-5] == '1'

    def FindEnd(self, Read):
        '''
        calculate the end of a mapping
        insertLength from start to end does not include the 5' or 3' soft clip part
        '''
        cigar, start = Read.split('\t')[5], int(Read.split('\t')[3])
        # ls = re.findall('\d+[MDIS]', cigar) # '4M1D72M' -> ['4M', '1D', '72M']
        # end = sum([int(re.findall('(\d+?)[MD]', item)[0]) for item in ls if 'S' not in item and 'I' not in item]) + start - 1
        end = sum([int(item) for item in re.findall('(\d+?)[MD]', cigar)]) + start - 1
        return start, end # int

    def AdjustCigar(self, cigar, insertLength):
        '''
        cigar: original cigar
        insertLength: 
        after trimming, new mapping end - new mapping start + 1, 
        not including the 5'end soft clip, equaling to sum(M|D) in new cigar

        adjust the cigar to the match the trimmed mapping
        e.g. original cigar is 1S5M1D1I6M, since the SEQ is shortened to 8nt, the returning new cigar is 1S5M1D1I1M
        Tips:
            len(new SEQ) = sum(M|I|S) in cigar, which should match the new cigar
        '''
        current_length = 0
        ls = []
        
        for item in zip(re.findall('\d+([MDIS])', cigar), re.findall('(\d+?)[MDIS]', cigar)):
            if item[0] not in ['I', 'S']:
                if current_length + int(item[1]) < insertLength:
                    ls.append('{}{}'.format(item[1], item[0]))
                    current_length = current_length + int(item[1])
                else:
                    ls.append('{}{}'.format(insertLength-current_length, item[0]))
                    break
            else:
                ls.append('{}{}'.format(item[1], item[0]))
        return ''.join(ls) # new cigar

    def AdjustMD(self, MD, insertLength):
        '''
        MD: original MD tag, e.g. MD:Z:95
        insertLength: 
        after trimming, new mapping end - new mapping start + 1, 
        not including the 5'end soft clip, equaling to sum(M|D) in new cigar

        adjust the MD tag to the match the trimmed mapping
        
        Tips:
            The soft clip part is counted for len(SEQ)=length but not shown in MD
            The insertion (I in cigar, base in read not in Reference) is not shown in MD
            For deletion (D in cigar, base in Reference not in read) e,g 2D in cigar, labeled as ^[ATCG][ATCG] in MD
        '''
        if insertLength <=0: # complete overlapping
            return 'MD:Z:0'
        else:
            MD = MD.split(':')[-1]
            current_length = 0
            ls = []
            for item in regex.findall('\d{1,}(?=[ATCGN])|[ATCGN]|\d{1,}|\^[ATCGN]+(?=\d)', MD):
                # '10T0T9A11G5T0G16T11^AG22' -> ['10', 'T', '0', 'T', '9', 'A', '11', 'G', '5', 'T', '0', 'G', '16', 'T', '11', '^AG', '22']
                if item.startswith('^'): # corresponding to D in cigar, bases in Reference but deleted in read
                    ls.append(item)
                    current_length = current_length + len(item)-1 # -1 since len contains '^'
                elif item.isdigit():
                    if int(item)<insertLength-current_length:
                        ls.append(item)
                        current_length = current_length+int(item)
                    else:
                        ls.append(insertLength-current_length)
                        break
                else:
                    if current_length<insertLength:
                        ls.append(item)
                        current_length +=1
                    else:
                        break
            return 'MD:Z:{}'.format(''.join([str(x) for x in ls])) # MD:Z:0C49

    def TrimOverlapping(self, Read1, Read2, count):
        '''
        Determine whether a pair of reads overlap
        count: number of overlapping read pair
        logic:
        If gtf_end of forward mapping read >= gtf_start of reverse mapping read, the pair of reads overlap.
            For overlapping read pair, trim the 3'end of the forward mapping read to remove the overlapping region;
                Trim SEQ, QUAL, adjust the cigar and MD tag accrodingly
                Leave the 5'end soft clip in the forward mapping.
        If not overlap, return Read1, Read2
        Cigar has other than 'MIDS' will not be trimmed.

        Note:
        For the complete overlapping, if leave the col10 and col11 as None, sam2bam conversion by samtools will add '*' accordingly,
        therefore even the sam file is not corret (having empty col10 and col11), the bam file format is correct.
        '''
        Read1_start, Read1_end = self.FindEnd(Read1)
        Read2_start, Read2_end = self.FindEnd(Read2)
        temp1, temp2 = Read1.split('\t'), Read2.split('\t')

        # To speed up, do not check [XHP] in cigar for now
        # if re.findall('[XHP=]', temp1[5]) or re.findall('[XHP=]', temp2[5]): # cigar not accepted
        #     print("Skip the analysis for this pair of read since Cigar containing string not MIDS: {}".format(temp1[0]))
        #     return Read1, Read2, count
        # else: # cigar has only 'MIDS'
        if not self.Flag16(temp1[1]): # Read1 forward mapping
            if Read1_end >= Read2_start: # overlapping
                Read1_end_new = Read2_start - 1

                # insertLength is the length of trimmed mapping part, not including the 5'end soft clip
                insertLength = Read1_end_new - Read1_start + 1

                if insertLength <=0: # complete overlapping: Read1_start<=Read2_start
                    Read1_cigar_new = '0M'
                    Read1_SEQ_new, Read1_QUAL_new = '*', '*'
                else: # not complete overlapping
                    # remove 3'end soft clip since trimming is from 3'end
                    if re.findall('\d+?S$', temp1[5]):
                        softclip_3_length = int(re.findall('(\d+?)S$', temp1[5])[0])
                        temp1[5] = temp1[5].replace(re.findall('\d+?S$', temp1[5])[0], '')
                        temp1[9] = temp1[9][:-1*softclip_3_length]
                        temp1[10] = temp1[10][:-1*softclip_3_length]
                    
                    # 5'end soft clip is included in the trimmed mapping
                    Read1_cigar_new = self.AdjustCigar(temp1[5], insertLength)
                    SEQlength = sum([int(item) for item in re.findall('(\d+)[MIS]', Read1_cigar_new)])
                    Read1_SEQ_new = temp1[9][:SEQlength]
                    Read1_QUAL_new = temp1[10][:SEQlength]
                
                Read1_new = temp1[:5]
                Read1_new.extend([Read1_cigar_new, temp1[6], temp1[7], temp1[8], Read1_SEQ_new, Read1_QUAL_new])
                
                # add tags
                for tag in temp1[11:]:
                    if tag.startswith('MD'):
                        tag = self.AdjustMD(tag, insertLength)
                    Read1_new.append(tag)
                
                # add custom XO tag to label the trimmed mapping
                Read1_new.append('XO:Z:readtrim')
                temp2.append('XO:Z:matetrim')

                count += 1
                return '\t'.join(Read1_new), '\t'.join(temp2), count

            else:
                return Read1, Read2, count

        else: # Read2 forward mapping
            if Read2_end >= Read1_start: # overlapping
                Read2_end_new = Read1_start - 1

                # insertLength is the length of trimmed insert, not including the 5'end soft clip
                insertLength = Read2_end_new - Read2_start + 1

                if insertLength <=0: # complete overlapping: Read2_start<=Read1_start
                    Read2_cigar_new = '0M'
                    Read2_SEQ_new, Read2_QUAL_new = '*', '*'
                else: # not complete overlapping
                    # remove 3'end soft clip since trimming is from 3'end
                    if re.findall('\d+?S$', temp2[5]):
                        softclip_3_length = int(re.findall('(\d+?)S$', temp2[5])[0])
                        temp2[5] = temp2[5].replace(re.findall('\d+?S$', temp2[5])[0], '')
                        temp2[9] = temp2[9][:-1*softclip_3_length]
                        temp2[10] = temp2[10][:-1*softclip_3_length]
                    
                    # 5'end soft clip is included in the trimmed mapping
                    Read2_cigar_new = self.AdjustCigar(temp2[5], insertLength)
                    SEQlength = sum([int(item) for item in re.findall('(\d+)[MIS]', Read2_cigar_new)])
                    Read2_SEQ_new = temp2[9][:SEQlength]
                    Read2_QUAL_new = temp2[10][:SEQlength]
                
                Read2_new = temp2[:5]
                Read2_new.extend([Read2_cigar_new, temp2[6], temp2[7], temp2[8], Read2_SEQ_new, Read2_QUAL_new])
                
                # add tags
                for tag in temp2[11:]:
                    if tag.startswith('MD'):
                        tag = self.AdjustMD(tag, insertLength)
                    Read2_new.append(tag)
                
                # add custom XO tag to label the trimmed mapping
                Read2_new.append('XO:Z:readtrim')
                temp1.append('XO:Z:matetrim')

                count += 1
                return '\t'.join(temp1), '\t'.join(Read2_new), count                    

            else:
                return Read1, Read2, count

    def GoThroughSam(self, input_sam, output_sam):
        '''
        Go through read pairs and trim the forward mapping read if overlap
        '''
        output = open(output_sam, 'w')
        with open(input_sam) as f:
            line = f.readline().strip()
            while line.startswith('@'):
                print(line, file=output)
                line = f.readline().strip()
            
            # here self.output_sam is the final output, not the temporary output_sam used here
            customHeader = ['@PG', 'ID:TrimOverlappingReadPair', 'VN:{}'.format(__version__.replace('.','')), \
                            'CL:{}'.format('python TrimOverlappingReadPair.py --input {} --output {}'.format(self.input_sam, self.output_sam))]
            print('\t'.join(customHeader), file=output)
            
            count = 0 # number of overlapping read pair
            q = deque(maxlen=2)
            q.append(line)
            line = f.readline().strip() # need strip() to remove the '\n' at the end
            q.append(line)
            while line:
                if q[0].split()[0] == q[1].split()[0]: # for reads in pair, the IDs are the same
                    Read1, Read2, count = self.TrimOverlapping(q[0], q[1], count)
                    print(Read1, file=output)
                    print(Read2, file=output)
                    line = f.readline().strip()
                    q.append(line)
                    line = f.readline().strip()
                    q.append(line)
                else: # Remove unpaired read from output
                    line = f.readline().strip()
                    q.append(line)
        output.close()
        print("{} overlapping Read pairs are trimmed.".format(count))
        return output_sam
    
    def RemoveEmptyMapping(self, input_sam, output_sam):
        '''remove the entries having no matching in cigar, which will raise error for bedtools'''
        output = open(output_sam, 'w')
        with open(input_sam) as f:
            line = f.readline().strip()
            while line.startswith('@'):
                print(line, file=output)
                line = f.readline().strip()
            if line.split()[9] != '*':
                print(line, file=output)
            for line in f:
                if line.split()[9] != '*':
                    print(line.strip(), file=output)
        output.close()
        return output_sam

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input sam file, sorted by name', dest='input_sam', action='store')
    parser.add_argument('--output', help='output sam file', dest='output_sam', action='store')
    parser.add_argument('--removeEmptyMapping', help='Add this option to remove empty mapping', action='store_true', \
                        dest='Flag', required=False)
                      
    args = parser.parse_args()

    if DetectSort(args.input_sam):
        TrimOverlapping(args.input_sam, args.output_sam, args.Flag)
    else:
        quit()