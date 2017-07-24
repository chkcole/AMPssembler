#!/usr/bin/python
"""
Created on Thu Jun 11 13:22:50 2015

@author: charles
"""
import sys
import argparse
import math
import collections
import string
out_2_handle = open("Failed_Amplicons.txt","w")
def main(args):
    file_handle = None
    try:
        file_handle = open(sys.argv[1],"rU")
    except:
        sys.stderr.write("Read file was not found\n")
        exit(1)
    for read_collection,Current_ID in Process_Reads(file_handle,33,150):
        generate_contigs(read_collection,Current_ID,15)
        

class contig(object):
    """
    This class described an object which will hold a contiguous piece of DNA 
    as well as the quality scores and counts for each base in the 5' to 3' 
    direction. It includes methods for creating a new contig, extending the contig
    and printing the contig object.
    """
    def __init__(self,R = "",Q = list() ,C = list()):
        self.R = R#Nucleotide sequence
        self.Q = Q#Phred quality score of each base
        self.C = C#Counts for each base
    def extend(self,R,Q,C):
        self.R += R
        self.Q.append(Q)
        self.C.append(C)
    def __repr__(self):
        return(self.R+"\n"+" ".join(map(str,self.Q))+"\n"+" ".join(map(str,self.C)))
class Kmer(object):
    """
    This object represents metainformation associated with a kmer.
    This information includes the additive quality of the final base of the kmer,
    the counts of that kmer and a flag indicating if this kmer is found at the beginning
    or the end of the amplicon.
    """
    def __init__(self,Q = 0,C = 0,is_anchor = False):
        self.Q = Q#Quality Score of final base
        self.C = C#Counter
        self.A = is_anchor#flag indicating if the kmer occurs at the beginning or end of the amplicon
    def __iadd__(self, other):
        self.Q += other.Q
        self.C += other.C
        return self
    def __repr__(self):
        return(str(self.C)+"\t"+str(self.Q))

def generate_contigs(read_collection,UID,kmer_size = 15):
    """
    This function generates contigs from a collection of reads given a kmer size to use for kmer extension
    The read_collecion is an array of readPair objects, the UID is simply the name you want to give the contig when it's printed to stdout in fastq format
    the kmer_size determines the size of kmers used to extend the anchors
    """
    nuc_trans = string.maketrans("ATGC","TACG")
    F_hash = collections.defaultdict(Kmer)#stores kmers generated directly from the reads in 5' to 3' direction
    R_hash = collections.defaultdict(Kmer)#stores kmer generected by the reverse complement of the reads in 5' to 3' direction
    C_anchor_array = [] #This array stores reads anchored to the constant region of the IG transcript
    V_anchor_array = [] #This array stores reads anchored to the variable region of the IG transcript
    C_anchor = None; # This variable will hold a contg object defining the the end of the amplicon corresponding to the constant region
    V_anchor = None; # This variable will hold a contg object defining the the end of the amplicon corresponding to the variable region
    C_count = 0 # number of reads anchored to the constant region
    V_count = 0 # number of reads anchored to the variable region
    V_C_count = 0 #number of reads anchored to both the variable and constant region.
    for read_pair in read_collection: #This loops breaks down the reads into forward and reverse kmers and counts the number of reads anchored to different ends.
        for kmer_tuple in generate_kmers(read_pair.Seq_1.R,read_pair.Seq_1.Q,kmer_size):
            F_hash[kmer_tuple[0]] += kmer_tuple[1]
        RC = read_pair.RC_R2()
        for kmer_tuple in generate_kmers(RC[0],RC[1],kmer_size):
            F_hash[kmer_tuple[0]] += kmer_tuple[1]
        for kmer_tuple in generate_kmers(read_pair.Seq_2.R,read_pair.Seq_2.Q,kmer_size):
            R_hash[kmer_tuple[0]] += kmer_tuple[1]
        RC = read_pair.RC_R1()
        for kmer_tuple in generate_kmers(RC[0],RC[1],kmer_size):
            R_hash[kmer_tuple[0]] += kmer_tuple[1]

        if read_pair.V == False:
            C_anchor_array.append((read_pair.Seq_2.R,read_pair.Seq_2.Q))
            C_count += 1
        elif read_pair.C == False:
            V_anchor_array.append((read_pair.Seq_1.R,read_pair.Seq_1.Q))
            V_count += 1
        else:
            C_anchor_array.append((read_pair.Seq_2.R,read_pair.Seq_2.Q))
            V_anchor_array.append((read_pair.Seq_1.R,read_pair.Seq_1.Q))
            V_C_count += 1
            
    C_anchor = concensus(C_anchor_array)
    V_anchor = concensus(V_anchor_array)
    C_end_kmer = string.translate(C_anchor.R[:kmer_size][::-1],nuc_trans)
    V_end_kmer = string.translate(V_anchor.R[:kmer_size][::-1],nuc_trans)    
    F_hash[C_end_kmer].A = True
    R_hash[V_end_kmer].A = True
    Extended_C,complete = extend_contig(C_anchor,R_hash,kmer_size)
    if complete:
        format_output(string.translate(Extended_C.R[::-1],nuc_trans),Extended_C.Q[::-1],V_C_count,C_count,V_count,UID)
        return True
    Extended_V,complete = extend_contig(V_anchor,F_hash,kmer_size)
    if complete:
        format_output(Extended_V.R,Extended_V.Q,V_C_count,C_count,V_count,UID)
        return True
    out_2_handle.write(str(UID)+"\t"+str(V_C_count)+"\t"+str(C_count)+"\t"+str(V_count)+"\n")
    return False

def extend_contig(begin_contig,kmer_hash,kmer_size,max_len = 1000):
    """
    This function attempts to extend a contig in the 5' to 3' direction using a hash of kmers.
    """
    extendable_contig = begin_contig
    iter_count = 0
    while True:
        iter_count += 1
        if iter_count > max_len:
            return(extendable_contig,False)
        current_kmer_counts = 0
        path_count = 0
        Qual= 0
        nuc = ""
        kmer_to_extend = extendable_contig.R[-(kmer_size - 1):]
        for n in "ATGCN":
            if kmer_hash[kmer_to_extend+n].Q > Qual:
                nuc = n
                path_count += 1
                Qual = kmer_hash[kmer_to_extend+n].Q
                current_kmer_counts = kmer_hash[kmer_to_extend+n].C
        if path_count  < 1:
            return(extendable_contig,False)
        else:
            extendable_contig.extend(nuc,Qual,current_kmer_counts)
            if kmer_hash[kmer_to_extend+nuc].A == True:
                return(extendable_contig,True)
                
def concensus(tuple_array):
    """
    This function creates a contig out of an array of sequence, quality-score pairs by finding the concensus of those sequences.
    It determines the concensus by adding up the quality scores for each base at each position and selecting the base with the highest
    total quality score.
    """
    counts = []
    sequence = contig("",[],[])
    for i in range(len(tuple_array[0][0])):
        counts.append([[0,0],[0,0],[0,0],[0,0],[0,0]])
    for i in tuple_array:
        for j,n in enumerate(i[0]):
            if n == "A":
                counts[j][0][0]+= 1
                counts[j][0][1] += i[1][j]
            elif n == "T":
                counts[j][1][0] += 1
                counts[j][1][1] += i[1][j]
            elif n == "G":
                counts[j][2][0] += 1
                counts[j][2][1] += i[1][j]
            elif n == "C":
                counts[j][3][0] += 1
                counts[j][3][1] += i[1][j]
            else:
                counts[j][4][0] += 1
                counts[j][4][1] += i[1][j]
    for i in counts:
        maximum = 0
        favored_nuc = "N"
        count_of = 0
        other_qual = 0
        for j,n in enumerate("ATGCN"):
            if i[j][1] > maximum:
                other_qual += maximum
                favored_nuc = n
                maximum = i[j][1]
                count_of = i[j][0]
            else:
                other_qual += i[j][1] 
        sequence.extend(favored_nuc,maximum - other_qual,count_of)
    return(sequence)


def generate_kmers(seq,qual,kmer_size = 15):
    """
    This function generates the information needed to create a kmer object
    from a string of characters, a quality score array and a kmer size specification.
    """
    for i in range(len(seq) - kmer_size + 1):
		yield(seq[i:i+kmer_size],Kmer(qual[i+kmer_size - 1],1,False))

class readPair(object):
    """
    This class represents a paired-end read. It contain the base sequence and quality scores for each end.
    """
    def __init__(self,Seq_1,Seq_2,Qual_1,Qual_2,C = True,V = True,encoding = 33):
        self.C = C # flag indicating if the readpair is anchored in the constant region
        self.V = V # flag indicating if the readpair is anchored in the variable region
        self.nuc_trans = string.maketrans("ATGC","TACG")
        length_check = len(Seq_1)
        #assert all( x == length_check for x in map(len,[Seq_2,Qual_1,Qual_2]))
        self.Seq_1 = fastqRead(Seq_1,Qual_1,encoding)#forward read
        self.Seq_2 = fastqRead(Seq_2,Qual_2,encoding)#reverse read
    def RC_R1(self): #return reverse complement of read 1
        return(string.translate(self.Seq_1.R[::-1],self.nuc_trans),self.Seq_1.Q[::-1])
    def RC_R2(self): #return reverse complement of read 2
        return(string.translate(self.Seq_2.R[::-1],self.nuc_trans),self.Seq_2.Q[::-1])
    def __repr__(self):
        return(self.Seq_1.R+"\t"+self.Seq_2.R+"\t"+str(self.C)+"\t"+str(self.V))

class fastqRead(object):
    """
    This object represents a single-ended read. It contains the base sequence and quality score
    """
    def __init__(self,seq,qual,Phred= 33):
        self.R = seq #nucleotide seuquence
        if Phred == 33:
            self.Q = map(lambda x: x - 33, map(ord,qual)) #quality score
            assert all(x >= 0 for x in self.Q)
        elif Phred == 64:
            self.Q = map(lambda x: x - 64, map(ord,qual)) #quality score
            assert all(x >= 0 for x in self.Q)
        elif Phred == "log":
            self.Q = qual #quality score
    def __repr__(self):
        return(self.R)

def Process_Reads(file_handle, encoding = 33,length = None):
    """
    This function reads a file containing paired end reads as well as anchoring and
    barcode information and uses them to generate readpair information.
    """
    Current_ID = ""
    Current_Collection = []
    line = ""
    for line in file_handle:
        try:
            la = line.split("\t")
            C = True
            V = True
            if la[2] == "":
                V = False
            elif la[3] == "":
                C = False
            
            if Current_ID != la[1]:
                if len(Current_Collection) > 0:
                    yield(Current_Collection,Current_ID)
                    Current_ID = la[1]
                    Current_Collection =[generateReadPairFromFile(la[4],la[5],C,V, encoding,length)]
                else:
                    Current_ID = la[1]
                    Current_Collection =[generateReadPairFromFile(la[4],la[5],C,V,encoding,length)]
            else:
                Current_Collection.append(generateReadPairFromFile(la[4],la[5],C,V, encoding,length))
        except:
            sys.stderr.write("Bad Read")
    if len(Current_Collection) > 0:
        yield(Current_Collection,Current_ID)                       
def generateReadPairFromFile(line_1,line_2,C,V,encoding = 33,length = 150):
    np = line_1.split("~~~~")
    qp = line_2.split("~~~~")
    if length == None:
        length = len(np[0])
    if length > len(np[0]):
        length = len(np[0])
    return(readPair(np[0][:length],np[1][:length],qp[0][:length],qp[1][:length],C,V,encoding))

def format_output(seq,qual,both,C,V,UID):
    print("@"+str(UID)+"_#"+str(both)+"_#"+str(C)+"_#"+str(V))
    print(seq)
    print("+")
    print("".join(map(chr,map(lambda x: max(x,2),map(lambda x: x + 33,map(lambda x: min(x,90),qual))))))
if(__name__ == "__main__"):
    sys.exit(main(sys.argv))