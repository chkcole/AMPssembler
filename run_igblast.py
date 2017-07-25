#!/usr/bin/python
import sys
import os

# Set IgBLAST parameters

path=sys.argv[1]
infile_name=path+'/allfiltered_numbered.fasta'
outfile_name=path+'/igblast_out.txt'

germline_db_V = '/home/charles/Documents/Ig_Blast/Vollmers_setup/Igblast/ncbi-igblast-1.4.0/germline/IMGT_V_formated'
germline_db_D = '/home/charles/Documents/Ig_Blast/Vollmers_setup/Igblast/ncbi-igblast-1.4.0/germline/IMGT_D_formated'
germline_db_J = '/home/charles/Documents/Ig_Blast/Vollmers_setup/Igblast/ncbi-igblast-1.4.0/germline/IMGT_J_formated'

igblast_bin_path = '/home/charles/Documents/Ig_Blast/ncbi-igblast-1.4.0-src/c++/ReleaseMT/bin/igblastn'
igblast_dir = '/home/charles/Documents/Ig_Blast/Vollmers_setup/Igblast/'	# must contain internal data directory

cwd = os.getcwd() # Get current working directory
os.chdir(igblast_dir)

#cmd = igblast_bin_path + ' -germline_db_V '+ germline_db_V+ ' -germline_db_D '+ germline_db_D+ ' -germline_db_J ' + germline_db_J +  ' -query ' + infile_name + ' -domain_system imgt -outfmt 3' + ' -out ' + outfile_name + ' -num_alignments_J=1 -num_alignments_D=1 -num_alignments_V=1 -evalue=1e-200'

cmd = igblast_bin_path  + ' -germline_db_V '+ germline_db_V+ ' -germline_db_D '+ germline_db_D+ ' -germline_db_J ' + germline_db_J +  ' -query ' + infile_name + ' -domain_system imgt -outfmt 3' + ' -out ' + outfile_name + ' -num_alignments_J=1 -num_alignments_D=1 -num_alignments_V=1 -evalue=1e-200'

os.system(cmd)



