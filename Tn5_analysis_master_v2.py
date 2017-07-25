#! /usr/bin/python
import os
import sys

for folder in open(sys.argv[1]):
  path=folder.strip().split('\t')[0]+'/'
  Chain_Type=folder.strip().split('\t')[1]
  print path
  for subdirectory in ['A','B','Uncut']:
    print(subdirectory)
    infile1_list=[]
    infile2_list=[]
    subpath=path+'/'+subdirectory+'/'
    os.system("gunzip "+subpath+"*"+Chain_Type+"*")
    os.system('rm '+subpath+'/Combined_R1.fastq')
    os.system('rm '+subpath+'/Combined_R2.fastq')
    for file1 in sorted(os.listdir(subpath)):
        if 'fastq' in file1:
            if 'gz' not in file1:
                if 'R1' in file1:
                    hypothetical_R2 = file1.replace('R1','R2')
                    if hypothetical_R2 in os.listdir(subpath):
                        infile1_list.append(subpath+file1)
                        infile2_list.append(subpath+hypothetical_R2)
                        print(subpath+file1,subpath+hypothetical_R2)
#    if len(infile1_list)==1:
#       print '1 file only'
#       infile1=infile1_list[0]
#       infile2=infile2_list[0]
#       print infile1, infile2

#    else:
    out_reads1=open(subpath+'Combined_R1.fastq','w')
    for element in infile1_list:
        print element
        for line in open(element):
           out_reads1.write(line)
    out_reads1.close()
    out_reads2=open(subpath+'Combined_R2.fastq','w')
    for element in infile2_list:
        print element
        for line in open(element):
          out_reads2.write(line)
    out_reads2.close()
    infile1=subpath+'Combined_R1.fastq'
    infile2=subpath+'Combined_R2.fastq'
    #os.system("gzip "+subpath+"*"+Chain_Type+"*fastq")
    os.system('java -jar trimmomatic-0.33.jar PE -phred33 '+ infile1 +' '+ infile2 +' '+subpath+'/R1_trimmed.txt.gz '+subpath+'/R1_trimmed_unpaired.txt.gz '+subpath+'/R2_trimmed.txt.gz '+subpath+'/R2_trimmed_unpaired.txt.gz ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10 MINLEN:20')
#    os.system('mv '+  subpath+'Combined_R1.fastq ' +subpath+'R1_trimmed.txt')
#    os.system('mv '+  subpath+'Combined_R2.fastq ' +subpath+'R2_trimmed.txt')
    #
    os.system('gunzip -f '+subpath+'/R*_trimmed.txt.gz')
#  out=open('summary_trimmomatic.txt','w')
  os.system('python index_reads_IGH_Tn5_Uncut.py '+path+' '+path+'/Uncut/ '+Chain_Type)
  os.system("gzip "+subpath+"*"+Chain_Type+"*")
  for cutoff_Uncut in ['5']:
      for cutoff_AB in ['2']:
          print 'cutoff Uncut ', cutoff_Uncut
          print 'cutoff AB ', cutoff_AB

          os.system('python index_reads_IGH_Tn5_Uncut_fuzzy.py '+path+' '+path+'/Uncut/ '+cutoff_Uncut + ' ' +Chain_Type)
          os.system('python index_reads_IGH_Tn5_AB_fuzzy.py '+path+' '+path+'/A/'+' '+path+'/B/ '+cutoff_AB )
          os.system('python cout_coverage.py '+path)
          print 'coverage calculation complete'
          os.system('python AMPssembler.py '+path+ '/indexed_reads_Uncut_AB_sorted.txt > '+path+'consensus_complete.txt')
          length=0
          for line in open(path+'consensus_complete.txt'):
               length+=1
          print length/4
          os.system('python combine_identicals_cluster.py '+ path)
          os.system("./run_igblast.py "+path)
  os.system("rm "+path+"/indexed*")
  os.system("rm "+path+"/A/*trimmed*")
  os.system("rm "+path+"/A/*Combined*")
  os.system("rm "+path+"/B/*trimmed*")
  os.system("rm "+path+"/B/*Combined*")
  os.system("rm "+path+"/Uncut/*trimmed*")
  os.system("rm "+path+"/Uncut/*Combined*")
#          out.write(str(cutoff_Uncut)+'\t'+str(cutoff_AB)+'\t'+str(length/4)+'\n')


