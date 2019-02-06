import re,os,sys,subprocess

#tissue_list = ["breast","endometrium","lung","ovary","pancreas","stomach","urinary_tract"]
tissue_list = ["breast"]
path = '/mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/output_ccle/ISO_module/'
IsoExon = '/mnt/isilon/xing_lab/shens/ISO_ReRun/CCLE/output_ccle/ISO_module/G20460.COR-L24.2.bam.sorted.IsoExon'

#CCLE_EMonly.txt
ifile=open('CCLE_EMonly.txt');
ilines=ifile.readlines();
bam=[];tissue=[];em=[];
for i in ilines:
	elements=re.findall('[^\t\n]+',i)
	#print(elements);
	bam.append(elements[1]);
	tissue.append(elements[9]);
	em.append(elements[13]);

#copy each file into the folder for each tissue
for i in tissue_list:
	subprocess.call("mkdir "+i, shell=True);
	subprocess.call("mkdir "+i+'/ISO_module', shell=True);
	subprocess.call("mkdir "+i+'/EM_out', shell=True);
	#sample list files
	ofile=open(i+'/input.list','w');
	flag = 0;e_count=0;m_count=0;
	e_list='';m_list='';
	for j in range(len(bam)):
	
		if tissue[j]==i:
			#copy isomatrix
			if (flag == 0) & (em[j]=='E'):
				subprocess.call("cp "+IsoExon+' '+i+'/ISO_module/'+bam[j]+'.sorted.IsoExon', shell=True);
				flag = 1
			subprocess.call("cp "+path+bam[j]+'.sorted.IsoMatrix'+' '+i+'/ISO_module/', shell=True);
			#write EM status
			if em[j]=='E':
				e_count+=1;
				e_list+=bam[j]+'.sorted\n'
			if em[j]=='M':
				m_count+=1;
				m_list+=bam[j]+'.sorted\n'
	
	ofile.write('2\n'+str(e_count)+'\n'+e_list+str(m_count)+'\n'+m_list);
	ofile.close();
