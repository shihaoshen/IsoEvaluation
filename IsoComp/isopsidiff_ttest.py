#02/12/19, use the GLMM logit regression results as a gold standard
#02/13/19, instead of using the major isoform (isoform with the largest sum of psi), use the isoform with largest difference between replicate average
#not to use rpy2 (loading issue on chop)
#use conda activate turbo

#import re,os,sys,numpy,rpy2
import re,os,sys,numpy,subprocess

# #rpy initialize
# import rpy2.robjects as robjects
# from rpy2.robjects import IntVector
# from rpy2.robjects.packages import importr
# # import R's "base" package
# base = importr('base')
# # import R's "utils" package
# utils = importr('utils')
# lme4 = importr("lme4")

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

def ttest(input_inc, input_tot, input_group):
	ofile2=open('glmer_input.txt','w');
	ofile2.write('inc\ttot\tgroup\n');
	for i in range(len(input_inc)):
		ofile2.write(str(int(input_inc[i]))+'\t'+str(int(input_tot[i]))+'\t'+str(input_group[i])+'\n')
	ofile2.close();
	subprocess.call("Rscript ttest.R ", shell=True)
	ifile2=open('ttest_output.txt');
	ifile2.readline();iline=ifile2.readline();
	elements=re.findall('[^ \t\n]+',iline);
	ifile2.close();
	print(elements);
	return(elements[1]);
	
ofile=open(sys.argv[1],'w');
ofile.write('ID\tFluxPsiDiff\tNumIsoform\tType\tRegP\tFluxPsiG1\tFluxPsiG2\tTotalCountG1\tTotalCountG2\n');
iline = 1;
while iline:
	psi1=[[],[],[]]; psi2=[[],[],[]];
	count1=[0,0,0]; count2=[0,0,0];
	psi1avg=0; psi2avg=0;
	count1avg=0; count2avg=0;
	p = 1;
	#for each isoform, record the count in three replicates
	count1_iso=[];count2_iso=[];
	
	#get the average psi value for the first group, replicate i
	for i in range(len(list1)):
		iline = ifile1[i].readline();
		elements=re.findall('[^ \t\n]+',iline);
		if len(elements)==0:
			continue;
		
		#get the number of isoforms
		num_isoform = int(elements[1]);
		
		#initialize the count record for each isoform
		if i == 0:
			for j in range(num_isoform):
				count1_iso.append([]); count2_iso.append([]);
		asm = elements[0];
		for j in range(num_isoform):
			psi1[i].append(float(elements[j+2]));
		#get the count of each isoform
		for j in range(num_isoform):
			count1[i]+=float(elements[j+num_isoform+2]);
			#for each isoform j, fill in the replicate count for replicate i
			count1_iso[j].append(float(elements[j+num_isoform+2]));
			
	for i in range(len(list1)):
		if i<len(psi1):
			psi1avg+=numpy.array(psi1[i])
			count1avg+=numpy.array(count1[i])
	
	#get the average psi value for the second group
	for i in range(len(list2)):
		iline = ifile2[i].readline();
		elements=re.findall('[^ \t\n]+',iline);
		if len(elements)==0:
			continue;
		asm = elements[0];
		for j in range(num_isoform):
			psi2[i].append(float(elements[j+2]));
		#get the count of each isoform
		for j in range(num_isoform):
			count2[i]+=float(elements[j+num_isoform+2]);
			#for each isoform j, fill in the replicate count for replicate i
			count2_iso[j].append(float(elements[j+num_isoform+2]));
			
	for i in range(len(list2)):
		if i<len(psi2):
			psi2avg+=numpy.array(psi2[i])
			count2avg+=numpy.array(count2[i])
	
	if len(elements)>0:
		psi1avg=psi1avg/3;psi2avg=psi2avg/3;
		#get the major isoform 
		#the majority isoform with the largest total psi in both sample groups	
		psisum = psi1avg + psi2avg;
		psiindex = numpy.argmax(psisum);
		psi1major = psi1avg[psiindex]; psi2major = psi2avg[psiindex];
		
		count1avg=count1avg/3;count2avg=count2avg/3;
		#get the abs difference
		psiavg = abs(psi1avg - psi2avg);
		#get the isoform with the largest difference between replicate average
		psiindex = numpy.argmax(psiavg);
		psi1diff = psi1avg[psiindex]; psi2diff = psi2avg[psiindex];
		#get the count from the isoform with the largest difference
		#for each replicate, extract the count for the isoform with largest difference
		input_inc = []; input_tot = []; input_group = [];
		for i in range(len(list1)):
			input_inc.append(count1_iso[psiindex][i]);
			input_tot.append(count1[i]);
			input_group.append(0);
		for i in range(len(list1)):
			input_inc.append(count2_iso[psiindex][i]);
			input_tot.append(count2[i]);
			input_group.append(1);
		#print('psiavg');print(psiavg);print(psi1avg);print(psi2avg);print(input_inc);print(input_tot);			
			
		# remove rpy2 calling
		# rpyinc = IntVector(input_inc); rpytot = IntVector(input_tot);
		# group = base.gl(2, len(list1), labels=["Ctl","Trt"])
		# res = lme4.glmer("cbind(rpyinc, rpytot - rpyinc) ~ group + (1 | group), family = binomial")
		# res2 = base.summary(res)
		# print(res2.rx2('coefficients'))
		
		#get the count for each isoform
		countavg = (count1avg + count2avg)/2
		psimax = [];
		#print('psiavg');print(psi1avg);print(psi2avg);
		if countavg >= count_cutoff:
			psimax=str(max(psiavg))
			#p value from logistic regression of the isoform with largest difference
			if len(sys.argv)>2:
				if sys.argv[2]=='1':
					p = ttest(input_inc, input_tot, input_group);
		else:
			psimax='NA'; psi1major='NA'; psi2major='NA'; p='NA';
		#print(psiavg);print('psidiff');print('\n');
		ofile.write(asm+'\t'+str(psimax)+'\t'+str(len(psiavg))+'\t'+elements[-1]+'\t'+str(p)+'\t'+str(psi1major)+'\t'+str(psi2major)+'\t'+str(count1avg)+'\t'+str(count2avg)+'\n');
ofile.close()
