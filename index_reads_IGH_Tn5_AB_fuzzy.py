from __future__ import division
import time
import sys
import os
import distance

path_root=sys.argv[1]+'//'
path1=sys.argv[2]+'//'
path2=sys.argv[3]+'//'
cutoff=int(sys.argv[4])

def find_match(Library,index, split_index1,split_index2,split_index3):

        index1=index[:6]
        index2=index[6:12]
        index3=index[12:18]
        
        query_set=set()
        for subindex,split_index in [(index1,split_index1),(index2,split_index2),(index3,split_index3)]:
             try:
                 for element in split_index[subindex]:
                     query_set.add(element)
             except:
                 pass

        min1=cutoff
        index_number=''
        for entry in list(query_set):
            dis=distance.hamming(index,entry[0])
            if dis<=min1:
                     index_number=entry[1]
                     min1=dis
  
        if index_number!='':                    
            out=open(path_root+'indexed_reads_Uncut_AB.txt','a') 
            if Library=='A':
                out.write(Library+'\t'+str(index_number)+'\t\t'+index+'\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
            if Library=='B':
                out.write(Library+'\t'+str(index_number)+'\t'+index+'\t\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
            out.close()

        else:
            out=open(path_root+'indexed_reads_Orphan_AB.txt','a') 
            if Library=='A':
                out.write(Library+'\t'+str(index_number)+'\t'+index+'\t\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
            if Library=='B':
                out.write(Library+'\t'+str(index_number)+'\t\t'+index+'\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
            out.close()




def gather_indices(path_root):
    infile=path_root+'/indexed_reads_Uncut_mismatch_sorted.txt'
    split_index1={}
    split_index2={}
    split_index3={}
    split_index4={}
    split_index5={}
    split_index6={}

    index_dict={}
    molecule_dict={}
    for line in open(infile):
        out=open(path_root+'indexed_reads_Uncut_AB.txt','a') 
        out.write(line)
        out.close()
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

        for subindex,split_index in [(index1,split_index1),(index2,split_index2),(index3,split_index3)]:
            try:
                split_index[subindex].append((relevant_index[:18],index_number))
            except:
                split_index[subindex]=[]
                split_index[subindex].append((relevant_index[:18],index_number))

        for subindex,split_index in [(index4,split_index4),(index5,split_index5),(index6,split_index6)]:
            try:
                split_index[subindex].append((relevant_index[18:],index_number))
            except:
                split_index[subindex]=[]
                split_index[subindex].append((relevant_index[18:],index_number))
                

    return split_index1,split_index2,split_index3,split_index4,split_index5,split_index6

def file_len(fileName):
        i=0
        for line in open(fileName):
            i+=1
            
        return i



out=open(path_root+'indexed_reads_Uncut_AB.txt','w')
out.close()

length1=file_len(path1+'/R1_trimmed.txt')
infile1_1=path1+'/R1_trimmed.txt'
infile1_2=path1+'/R2_trimmed.txt'
length2=file_len(path2+'/R1_trimmed.txt')
infile2_1=path2+'/R1_trimmed.txt'
infile2_2=path2+'/R2_trimmed.txt'

split_index1,split_index2,split_index3,split_index4,split_index5,split_index6=gather_indices(path_root)



in1=open(infile1_1,'r')
in2=open(infile1_2,'r')
start=time.time()




x=4
A_total=0
while x<=length1:     
        
        a=in1.readline()
        b=in1.readline()
        c=in1.readline()
        d=in1.readline()
    
        e=in2.readline()
        f=in2.readline()
        g=in2.readline()
        h=in2.readline()

        

        position=a.strip().split(' ')[0]

        index2=f[0:18]
        read1=b.strip()+'~~~~'+f[18:].strip()
        qual1=d.strip()+'~~~~'+h[18:].strip()

       
        find_match('A',index2,split_index4,split_index5,split_index6)

        x+=4


in1=open(infile2_1,'r')
in2=open(infile2_2,'r')
x=4
B_total=0
while x<=length2:     

        a=in1.readline()
        b=in1.readline()
        c=in1.readline()
        d=in1.readline()
    
        e=in2.readline()
        f=in2.readline()
        g=in2.readline()
        h=in2.readline()

        position=a.strip().split(' ')[0]

        index1=b[0:18]
        read1=b[18:].strip()+'~~~~'+f.strip()
        qual1=d[18:].strip()+'~~~~'+h.strip()
    
  
        find_match('B',index1,split_index1,split_index2,split_index3)

        x+=4

print A_total,B_total
os.system('sort -k2,2n -k7,7n '+path_root+'/indexed_reads_Uncut_AB.txt >'+path_root+'//indexed_reads_Uncut_AB_sorted.txt')     

       
    



    
