import sys
path=sys.argv[1]
target=path+'//consensus_complete.txt'


out=open(path+'/allfiltered.fasta','w')
out1=open(path+'/allfiltered_numbered.fasta','w')
out_coverage=open(path+'/coverage.txt','w')

length=0
for line in open(target):
    length+=1
x=0
a=open(target)
reads={}
reads_numbers={}

coverage_dict={}
for x in xrange(0,20000,1):
    coverage_dict[x]=0


x=0
print length
while x<length:
    b=a.readline()
    c=a.readline()
    k=a.readline()
    l=a.readline()

 #   counts=int(b.split('_')[8])
 #   coverage_dict[counts]+=1


    if c in reads:
        reads[c]=reads[c]+b.strip().split('@')[1]+','
        reads_numbers[c]=reads_numbers[c]+1

    else:
        reads[c]=b.strip().split('@')[1]+','
        reads_numbers[c]=1

    x+=4

print len(reads)
for sequence, names in reads.iteritems():
    out.write('>'+names+'\n')
    out.write(sequence)

x=1
for sequence, names in reads_numbers.iteritems():
    out1.write('>'+str(x)+'_'+str(names)+'\n')
    out1.write(sequence)
    x+=1

for x in xrange(0,2000,1):
    out_coverage.write(str(coverage_dict[x])+'\t')
