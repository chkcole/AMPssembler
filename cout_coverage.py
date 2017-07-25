from __future__ import division
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
path=sys.argv[1]

fig = plt.figure(figsize=(4,4),frameon=True)

def boxplot_pretty(abundance,position,color):
                 ax=plt.gca()
                 bp=plt.boxplot(abundance,notch=0, positions=[position],sym='',widths=1)

                 plt.setp(bp['boxes'], color='black',linewidth=0.2,alpha=1)
                 plt.setp(bp['whiskers'], color='black',linestyle='-', linewidth=0.5,alpha=1)
                 plt.setp(bp['fliers'], color='black', marker='.', linewidth=0)
                 plt.setp(bp['medians'], color='black', linewidth=0.5,alpha=1)
                 plt.setp(bp['caps'], color='black', linewidth=0.5,alpha=1)

                 # Now fill the boxes with desired colors
                 boxColors = [color]
                 box = bp['boxes'][0]
                 boxX = []
                 boxY = []

                 for j in numpy.arange(0,len(box.get_xdata()),1):
                      boxX.append(box.get_xdata()[j])
                 for j in numpy.arange(0,len(box.get_ydata()),1):
                      boxY.append(box.get_ydata()[j])
                 boxCoords = zip(boxX,boxY)
                 boxPolygon = Polygon(boxCoords, facecolor=boxColors[0],alpha=0.5,linewidth=0.5)
                 ax.add_patch(boxPolygon)
                 # Now draw the median lines back over what we just filled in


def plot_data(dict1,axis1):
    
    for element in dict1:
        list1=dict1[element]

        plt.axes(axis1,frameon=True)
        if len(list1)>3:
            boxplot_pretty(list1,element,'black')
            plt.ylim(1,100) 
            plt.xlim(0,100)

def extract_data(list1,axis1,plot):
    
    bins=xrange(0,1000,1)

    Uncut=[]
    A=[]
    B=[]
    for element in list1:
        Uncut.append(element[0])
        A.append(element[1])
        B.append(element[2])
    Uncut_Coverage=(numpy.average(Uncut), numpy.std(Uncut))    
    A_Coverage=(numpy.average(A), numpy.std(A))      
    B_Coverage=(numpy.average(B), numpy.std(B))  

    density,bins=numpy.histogram(Uncut,bins)
    if plot==True:
        plt.axes(axis1,frameon=True)
        plt.plot(bins[:-1],density, color='grey')
        plt.yscale('log')
        plt.ylim(1,100000) 
        plt.xlim(0,100)



    return Uncut_Coverage, A_Coverage,B_Coverage







dict_count={}
for line in open(path+'indexed_reads_Uncut_AB_sorted.txt'):
    a=line.strip().split('\t')
    UID_number=a[1]
    type1=a[0]
    try:
        bla=dict_count[UID_number]
        dict_count[UID_number][type1]+=1
    except:
        dict_count[UID_number]={}
        dict_count[UID_number]['0']=0
        dict_count[UID_number]['B']=0
        dict_count[UID_number]['A']=0
        dict_count[UID_number][type1]+=1






All=[]
Uncut_Only=[]
Uncut_and_A=[]
Uncut_and_B=[]
Total_Uncut=[]

coverage_uncut_A={}
coverage_uncut_B={}

for x in xrange(0,10000,1):
    coverage_uncut_A[x]=[]
    coverage_uncut_B[x]=[]

Total_U=0
Total_A=0
Total_B=0


for number in dict_count:
    
    Total_U+=dict_count[number]['0']
    Total_A+=dict_count[number]['A']
    Total_B+=dict_count[number]['B']


    Total_Uncut.append((dict_count[number]['0'],dict_count[number]['A'],dict_count[number]['B']))

    coverage_uncut_A[dict_count[number]['0']].append(dict_count[number]['A'])
    coverage_uncut_B[dict_count[number]['0']].append(dict_count[number]['B'])

    if dict_count[number]['A']==0:
        if  dict_count[number]['B']==0:
            Uncut_Only.append((dict_count[number]['0'],dict_count[number]['A'],dict_count[number]['B']))
        else:
            Uncut_and_B.append((dict_count[number]['0'],dict_count[number]['A'],dict_count[number]['B']))
    else:
        if  dict_count[number]['B']==0:
            Uncut_and_A.append((dict_count[number]['0'],dict_count[number]['A'],dict_count[number]['B']))
        else:
            All.append((dict_count[number]['0'],dict_count[number]['A'],dict_count[number]['B']))





Uncut_Only_Uncut_Coverage, Uncut_Only_A_Coverage,Uncut_Only_B_Coverage=extract_data(Uncut_Only,[0,00.6,1,0.18],False)
Uncut_and_A_Uncut_Coverage,Uncut_and_A_A_Coverage,Uncut_and_A_B_Coverage=extract_data(Uncut_and_A,[0,0.2,1,0.18],False)
Uncut_and_B_Uncut_Coverage,Uncut_and_B_A_Coverage,Uncut_and_B_B_Coverage=extract_data(Uncut_and_B,[0,0,1,0.18],False)
All_Uncut_Coverage,All_A_Coverage,All_B_Coverage=extract_data(All,[0,0.4,1,0.18],False)
Total_Uncut_Uncut_Coverage,Total_Uncut_B_Coverage,Total_Uncut_B_Coverage=extract_data(Total_Uncut,[0,0.66,1,0.3],True)
plot_data(coverage_uncut_A,[0,0,1,0.3])
plot_data(coverage_uncut_B,[0,0.33,1,0.3])
 
print 'Uncut', Total_U
print 'A', Total_A
print 'B', Total_B



print 'Total_Uncut',len(Total_Uncut),Total_Uncut_Uncut_Coverage,Total_Uncut_B_Coverage,Total_Uncut_B_Coverage
print 'Uncut_Only',len(Uncut_Only), Uncut_Only_Uncut_Coverage, Uncut_Only_A_Coverage,Uncut_Only_B_Coverage
print 'All',len(All),All_Uncut_Coverage,All_A_Coverage,All_B_Coverage
print 'Uncut_and_A',len(Uncut_and_A),Uncut_and_A_Uncut_Coverage,Uncut_and_A_A_Coverage,Uncut_and_A_B_Coverage
print 'Uncut_and_B',len(Uncut_and_B),Uncut_and_B_Uncut_Coverage,Uncut_and_B_A_Coverage,Uncut_and_B_B_Coverage
            




plt.savefig(path+'read_coverage.pdf',dpi=1200,bbox_inches='tight')  
        














        
