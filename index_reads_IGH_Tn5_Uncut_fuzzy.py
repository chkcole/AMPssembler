from __future__ import division
import time
import sys
import os
import distance

path_root=sys.argv[1]+'//'
path1=sys.argv[2]+'//'
cutoff=int(sys.argv[3])
Chain_Type=sys.argv[4]

IGHC_dict={}
IGHC_dict['GGGAAGTAGTCCTTGACCAGGCAGCCC']='IGHG'
IGHC_dict['GGGGAAGAAGCCCTGGACCAGGCA']='IGHA'
IGHC_dict['AAGTAGCCCGTGGCCAGGCAGCCCAGA']='IGHE'
IGHC_dict['TGGGTGGTACCCAGTTATCAAGCATGCCAGGA']='IGHD'
IGHC_dict['GAAGGAAGTCCTGTGCGAGGCAGCCAA']='IGHM'


IGLC_dict={}
IGLC_dict['AGTGTGGCCTTGTTGGCTTGAAGCTCCTC']='IGLC'
IGLC_dict['AGTGTGGCCTTGTTGGCTTGGAGCTCCTC']='IGLC'

IGKC_dict={}
IGKC_dict['GAAGATGAAGACAGATGGTGCAGCCACAGTT']='IGKC'


Chain_dicts={}
Chain_dicts['IGH']=IGHC_dict
Chain_dicts['IGL']=IGLC_dict
Chain_dicts['IGK']=IGKC_dict


IGHVL_dict={}
IGHVL_dict['ATG']='Leader'


def gather_indices(path_root):

    

    infile=path_root+'/indexed_reads_Uncut.txt'
    split_index1={}
    split_index2={}
    split_index3={}
    split_index4={}
    split_index5={}
    split_index6={}

    index_dict={}
    molecule_dict={}
    for line in open(infile):
        a=line.strip().split('\t')

        read_id=a[-1]

        index_number=a[1]
        relevant_index=a[2]+a[3]
        index_dict[relevant_index]=index_number

        index1=relevant_index[:6]
        index2=relevant_index[6:12]
        index3=relevant_index[12:18]
        index4=relevant_index[18:24]
        index5=relevant_index[24:30]
        index6=relevant_index[30:36]

        for subindex,split_index in [(index1,split_index1),(index2,split_index2),(index3,split_index3),(index4,split_index4),(index5,split_index5),(index6,split_index6)]:
            try:
                split_index[subindex].append((relevant_index,index_number))
            except:
                split_index[subindex]=[]
                split_index[subindex].append((relevant_index,index_number))
                

    for iterator in xrange(0,2,1):
      changes=0
      print 'iteration', iterator
      for index in index_dict:
        index1=index[:6]
        index2=index[6:12]
        index3=index[12:18]
        index4=index[18:24]
        index5=index[24:30]
        index6=index[30:36]
        query_set=set()
        for subindex,split_index in [(index1,split_index1),(index2,split_index2),(index3,split_index3),(index4,split_index4),(index5,split_index5),(index6,split_index6)]:
             for element in split_index[subindex]:
                 query_set.add(element)
        for entry in list(query_set):
            dis=distance.hamming(index,entry[0])
            if dis<=cutoff:
                if index_dict[entry[0]]!=index_dict[index]:
                     index_dict[entry[0]]=index_dict[index]
                     changes+=1
      print changes

    return index_dict


def file_len(fileName):
        i=0
        for line in open(fileName):
            i+=1
            
        return i


out=open(path_root+'indexed_reads_Uncut_mismatch.txt','w')
out.close()

length1=file_len(path1+'/R1_trimmed.txt')
infile1_1=path1+'/R1_trimmed.txt'
infile1_2=path1+'/R2_trimmed.txt'

compare_indexes=gather_indices(path_root)

in1=open(infile1_1,'r')
in2=open(infile1_2,'r')
Chain_dict=Chain_dicts[Chain_Type]
x=4
while x<=length1:     
        match=0
        a=in1.readline()
        b=in1.readline()
        c=in1.readline()
        d=in1.readline()
    
        e=in2.readline()
        f=in2.readline()
        g=in2.readline()
        h=in2.readline()

        position=a.strip().split(' ')[0]

        read_id=a.split()[0]

        index1=b[0:18]
        index2=f[0:18]
        read1=b[18:].strip()+'~~~~'+f[18:].strip()
        qual1=d[18:].strip()+'~~~~'+h[18:].strip()
        index=index1+index2    

        
        try:
            bla=IGHVL_dict[b[18:21]] 
            go=0 
            for entry in Chain_dict:
                length=len(entry)
                         
                dis=distance.hamming(entry,f[18:18+length])
                if dis<=2:
                    go=1

            if go==1:
                index_number=compare_indexes[index]
                out=open(path_root+'indexed_reads_Uncut_mismatch.txt','a') 
                out.write('0\t'+str(index_number)+'\t'+index1+'\t'+index2+'\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
                out.close()        
        except:
            pass


        x+=4
os.system('sort -k2,2n -k7,7n '+path_root+'/indexed_reads_Uncut_mismatch.txt >'+path_root+'//indexed_reads_Uncut_mismatch_sorted.txt')     

       
    



    
