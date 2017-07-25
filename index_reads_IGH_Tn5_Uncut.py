from __future__ import division
import sys
import os
import distance

path_root=sys.argv[1]+'//'
path1=sys.argv[2]+'//'
Chain_Type=sys.argv[3]

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



def file_len(fileName):
        i=0
        for line in open(fileName):
            i+=1
            
        return i


out=open(path_root+'indexed_reads_Uncut.txt','w')
out.close()

length1=file_len(path1+'/R1_trimmed.txt')
infile1_1=path1+'/R1_trimmed.txt'
infile1_2=path1+'/R2_trimmed.txt'

compare_indexes={}
compare_indexes1={}
compare_indexes2={}
in1=open(infile1_1,'r')
in2=open(infile1_2,'r')


missedA=0
missedB=0
totalA=0
totalB=0
totalComplete=0
x=4
Chain_dict=Chain_dicts[Chain_Type]
while x<=length1:     
        totalComplete+=1
        a=in1.readline()
        b=in1.readline()  #Omit the [:150] if you do not intend to truncate your reads
        c=in1.readline()
        d=in1.readline()
    
        e=in2.readline()
        f=in2.readline()
        g=in2.readline()
        h=in2.readline()

        position=a.strip().split(' ')[0]

        index1=b[0:18]
        index2=f[0:18]
        read1=b[18:].strip()+'~~~~'+f[18:].strip()
        qual1=d[18:].strip()+'~~~~'+h[18:].strip()
        index=index1+'_'+index2     

        

        try:
            bla=IGHVL_dict[b[18:21]]
            go=0
            for entry in Chain_dict:
                length=len(entry)

                dis=distance.hamming(entry,f[18:18+length])
                if dis<=2:
                    go=1
            if go==1:                
              try:
                UID_number=compare_indexes[index]


              except:
                compare_indexes[index]=int(x/4)
                compare_indexes1[index1]=int(x/4)
                compare_indexes2[index2]=int(x/4)
                UID_number=int(x/4)

              out=open(path_root+'indexed_reads_Uncut.txt','a')  
              out.write('0\t'+str(UID_number)+'\t'+index1+'\t'+index2+'\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
              out.close()
            

        except:
            pass


        x+=4


os.system('sort -k2,2n -k7,7n '+path_root+'/indexed_reads_Uncut.txt >'+path_root+'//indexed_reads_Uncut_sorted.txt')     

       
    



    
