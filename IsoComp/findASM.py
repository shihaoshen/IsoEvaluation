#this script find the ASM matching with SE events

import re,os,sys
from decimal import Decimal

#read the SE events
ifile=open(sys.argv[1]);
filelist=ifile.read();
filelist=re.findall('[^,\t\n\r]+',filelist)
ifile.close();
ifile=open(sys.argv[2]);
plist=ifile.read();
plist=re.findall('[^,\t\n\r]+',plist)
ifile.close();
Inc={}; Skp={};
type_list=['SE','A3SS','A5SS','MXE']
for i in range(len(filelist)):
	#read the list of p values for events
	ifile=open(plist[i]);
	ifile.readline();
	ilines=ifile.readlines();
	id2p={};
	for j in ilines:
		elements=re.findall('[^\t\n]+',j);
		id=elements[0];
		if i <= 2:
			id2p[id]="{:.2E}".format(Decimal(float(elements[18])));
		if i == 3:
			id2p[id]="{:.2E}".format(Decimal(float(elements[20])));    
	ifile.close();
	#read the list of events
	ifile=open(filelist[i]);
	ifile.readline();
	ilines=ifile.readlines();
	for j in ilines:
		elements=re.findall('[^\t\n]+',j);
		id=elements[0];
		if i == 0:
			this_inc1=elements[8]+'_'+elements[5];
			this_inc2=elements[6]+'_'+elements[9];
			this_skp1=elements[8]+'_'+elements[9];
			this_skp2=elements[8]+'_'+elements[9];
		if i == 1:
			if elements[4] == '+':
				this_inc1=elements[10]+'_'+elements[5];
				this_inc2=elements[10]+'_'+elements[5];
				this_skp1=elements[10]+'_'+elements[7];
				this_skp2=elements[10]+'_'+elements[7];
			if elements[4] == '-':
				this_inc1=elements[6]+'_'+elements[9];
				this_inc2=elements[6]+'_'+elements[9];
				this_skp1=elements[8]+'_'+elements[9];
				this_skp2=elements[8]+'_'+elements[9];
		if i == 2:
			if elements[4] == '+':
				this_inc1=elements[8]+'_'+elements[9];
				this_inc2=elements[8]+'_'+elements[9];
				this_skp1=elements[6]+'_'+elements[9];
				this_skp2=elements[6]+'_'+elements[9];
			if elements[4] == '-':
				this_inc1=elements[10]+'_'+elements[5];
				this_inc2=elements[10]+'_'+elements[5];
				this_skp1=elements[10]+'_'+elements[7];
				this_skp2=elements[10]+'_'+elements[7];
		if i == 3:
			this_inc1=elements[10]+'_'+elements[5];
			this_inc2=elements[6]+'_'+elements[11];
			this_skp1=elements[10]+'_'+elements[7];
			this_skp2=elements[8]+'_'+elements[11];
		if id in id2p:
			if this_inc1 in Inc:
				Inc[this_inc1].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Inc[this_inc1]=[type_list[i]+'_'+id+'_'+id2p[id]];
			if this_inc2 in Inc:
				Inc[this_inc2].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Inc[this_inc2]=[type_list[i]+'_'+id+'_'+id2p[id]];
			if this_skp1 in Skp:
				Skp[this_skp1].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Skp[this_skp1]=[type_list[i]+'_'+id+'_'+id2p[id]];
			if this_skp2 in Skp:
				Skp[this_skp2].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Skp[this_skp2]=[type_list[i]+'_'+id+'_'+id2p[id]];
	ifile.close();
	
#read the junctions for each ASM
ifile=open(sys.argv[3]);
ilines=ifile.readlines();
ofile=open(sys.argv[4],'w');
ofile.write('ASMID\tMatch\tSE\tA3SS\tA5SS\tMXE\tSEFDR5\tA3SSFDR5\tA5SSFDR5\tMXEFDR5\tPmin\n')
IncEvent = []; SkpEvent=[]; thisASM = '';
#Signficant cutoff for SE, A3SS, A5SS, MXE
FDRCut = [0.0000631470588356,0.0000289591477134,0.000663722761454,0.000086155609442]
for i in ilines:
	elements=re.findall('[^\t\n\r]+',i);
	if 'ASM#' in elements[0]:
		#write information for the preivous ASM
		thisEvent = list(set(IncEvent) & set(SkpEvent));
		if thisASM != '':
			ofile.write(thisASM+'\t'); thisoutput = '';
			event = [0,0,0,0];
			for j in thisEvent:
				elements2 = re.findall('[^_]+',j);
				for k in range(len(type_list)):
					if elements2[0] == type_list[k]:
						event[k] += 1;
			eventoutput = '';
			for j in event:
				eventoutput += '\t'+str(j);										
			p = [0,0,0,0]; pmin = 2;
			for j in thisEvent:
				thisoutput += j+',';
				elements2 = re.findall('[^_]+',j);
				for k in range(len(type_list)):
					if elements2[0] == type_list[k]:
						pmin = min(pmin, float(elements2[2]))
						print(float(elements2[2]));
						print(min(pmin, float(elements2[2])));
						print('===');
						if float(elements2[2]) <= FDRCut[k]:
							p[k] += 1;
			poutput = '';
			for j in p:
				poutput += '\t'+str(j);							
			if len(thisEvent) > 0:
				thisoutput = thisoutput[:-1];
			else:
				thisoutput = 'NA';
			if pmin == 2:
				ofile.write(thisoutput+eventoutput+poutput+'\t'+str('NA')+'\n');
			else:
				ofile.write(thisoutput+eventoutput+poutput+'\t'+str(pmin)+'\n');
		#write information for the current ASM
		thisASM = elements[0];	iscoord = 1;
		IncEvent = []; SkpEvent=[];
		continue;
	#read coordinates
	if iscoord == 1:
		iscoord = 0; juc1 = []; juc2 = [];
		for j in elements:
			jucs = re.findall('[^,]+',j);
			juc1.append(str(int(jucs[0])-1));
			juc2.append(jucs[1]);
	#read junctions
	else:
		for j in range(len(elements)-2):
			thisjuc=juc2[int(elements[j+1])]+'_'+juc1[int(elements[j+2])];
			if thisjuc in Inc:
				IncEvent+=Inc[thisjuc];
			if thisjuc in Skp:
				SkpEvent+=Skp[thisjuc];
