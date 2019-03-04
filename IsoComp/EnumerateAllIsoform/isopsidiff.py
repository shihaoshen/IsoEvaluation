#this script get the max difference in isoform ratios
import re,os,sys,numpy

#cutoff for read counts
count_cutoff = 10;

#group1 psi
list1=['PC3E-1Aligned.sort.bam.psi','PC3E-2Aligned.sort.bam.psi','PC3E-3Aligned.sort.bam.psi'];
psi1=[]; ifile1=[];
for i in list1:
	ifile1.append(open(i));
	psi1.append([]);

#group2	psi
list2=['GS689.LI-1Aligned.sort.bam.psi','GS689.LI-2Aligned.sort.bam.psi','GS689.LI-3Aligned.sort.bam.psi'];
psi2=[]; ifile2=[];
for i in list2:
	ifile2.append(open(i));
	psi2.append([]);

ofile=open(sys.argv[1],'w');
ofile.write('ID\tFluxPsiDiff\tNumIsoform\tType\tFluxPsiG1\tFluxPsiG2\n');
iline = 1;
while iline:
	psi1=[[],[],[]]; psi2=[[],[],[]];
	count1=[0,0,0]; count2=[0,0,0];
	psi1avg=0; psi2avg=0;
	count1avg=0; count2avg=0;
	#get the average psi value for the first group
	for i in range(len(list1)):
		iline = ifile1[i].readline();
		elements=re.findall('[^ \t\n]+',iline);
		if len(elements)==0:
			continue;
		asm = elements[0];
		for j in range(2,len(elements)-(len(elements) - 2 - 1)/2 - 1):
			psi1[i].append(float(elements[j]));
		for j in range(len(elements)-(len(elements) - 2 - 1)/2 + 1 - 1, len(elements) - 1):
			count1[i]+=float(elements[j]);
	for i in range(len(list1)):
		if i<len(psi1):
			psi1avg+=numpy.array(psi1[i])
			count1avg+=numpy.array(count1[i])
	#print(psi1avg);print('psi1avg');
	#get the average psi value for the second group
	for i in range(len(list2)):
		iline = ifile2[i].readline();
		elements=re.findall('[^ \t\n]+',iline);
		if len(elements)==0:
			continue;
		asm = elements[0];
		for j in range(2,len(elements)-(len(elements) - 2 - 1)/2 - 1):
			psi2[i].append(float(elements[j]));
		for j in range(len(elements)-(len(elements) - 2 - 1)/2 + 1 - 1, len(elements) - 1):
			count2[i]+=float(elements[j]);
	for i in range(len(list2)):
		if i<len(psi2):
			psi2avg+=numpy.array(psi2[i])
			count2avg+=numpy.array(count2[i])
	#print(psi2avg);print('psi2avg');
	if len(elements)>0:
		psi1avg=psi1avg/3;psi2avg=psi2avg/3;
		#get the major isoform 
		
		#the majority isoform with the largest total psi in both sample groups	
		psisum = psi1avg + psi2avg;
		psiindex = numpy.argmax(psisum);
		psi1major = psi1avg[psiindex]; psi2major = psi2avg[psiindex];
		
		count1avg=count1avg/3;count2avg=count2avg/3;
		#get the abs difference
		psiavg = abs(psi1avg - psi2avg)
		#get the count for each isoform
		countavg = (count1avg + count2avg)/2
		psimax = [];
		#print('psiavg');print(psi1avg);print(psi2avg);
		if countavg >= count_cutoff:
			psimax=str(max(psiavg))
		else:
			psimax='NA'; psi1major='NA'; psi2major='NA'
		#print(psiavg);print('psidiff');print('\n');
		ofile.write(asm+'\t'+str(psimax)+'\t'+str(len(psiavg))+'\t'+elements[-1]+'\t'+str(psi1major)+'\t'+str(psi2major)+'\n');
ofile.close()
