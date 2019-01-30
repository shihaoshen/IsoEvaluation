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
		id=re.findall('clu_\d+_NA',elements[4])[0];
		#print('');print(id);print(elements[3]);
		if 'NA' in elements[3]:
			id2p[id]="NA";
		else:
			id2p[id]="{:.2E}".format(Decimal(float(elements[3])));
	ifile.close();
  
for i in range(len(filelist)):
	#read the list of events
	ifile=open(filelist[i]);
	ifile.readline();
	ilines=ifile.readlines();
	for j in ilines:
		elements=re.findall('[^\t\n]+',j);
		more=re.findall('[^:]+',elements[0]);
		id = more[3];
		this_inc1=more[1];
		this_skp1=more[2];
		if id in id2p:
			if this_inc1 in Inc:
				Inc[this_inc1].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Inc[this_inc1]=[type_list[i]+'_'+id+'_'+id2p[id]];
			if this_skp1 in Skp:
				Skp[this_skp1].append(type_list[i]+'_'+id+'_'+id2p[id]);
			else:
				Skp[this_skp1]=[type_list[i]+'_'+id+'_'+id2p[id]];
	ifile.close();

print(Inc)
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
		#thisEvent = list(set(IncEvent) & set(SkpEvent));
		thisEvent = list(set(IncEvent) | set(SkpEvent));
		#print(thisEvent)
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
						pmin = min(pmin, float(elements2[4]))
						print(float(elements2[2]));
						print(min(pmin, float(elements2[4])));
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
			#thisjuc=juc2[int(elements[j+1])]+'_'+juc1[int(elements[j+2])];
			thisjuc2=juc2[int(elements[j+1])];
			thisjuc1=juc1[int(elements[j+2])];
			#print(thisjuc1);
			if thisjuc1 in Inc:
				IncEvent+=Inc[thisjuc1];
			if thisjuc1 in Skp:
				SkpEvent+=Skp[thisjuc1];
			if thisjuc2 in Inc:
				IncEvent+=Inc[thisjuc2];
			if thisjuc2 in Skp:
				SkpEvent+=Skp[thisjuc2];
