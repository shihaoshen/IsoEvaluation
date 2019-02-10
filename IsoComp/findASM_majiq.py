#this script find the ASM matching with SE events

import re,os,sys
from decimal import Decimal

#read the SE events
ifile=open(sys.argv[1]);
filelist=ifile.read();
filelist=re.findall('[^,\t\n\r]+',filelist)
ifile.close();
# ifile=open(sys.argv[2]);
# plist=ifile.read();
# plist=re.findall('[^,\t\n\r]+',plist)
# ifile.close();
Inc={}; Skp={};
type_list=['MJ']
for i in range(len(filelist)):
	id2p={};
	#read the list of p values for events
	# ifile=open(plist[i]);
	# ifile.readline();
	# ilines=ifile.readlines();
	# id2p={};
	# for j in ilines:
		# elements=re.findall('[^\t\n]+',j);
		# id=elements[0];
		# id2p[id]="{:.2E}".format(Decimal(float(elements[17])));
	# ifile.close();
	#read the list of events
	ifile=open(filelist[i]);
	ifile.readline();
	ilines=ifile.readlines();
	for j in ilines:
		elements=re.findall('[^\t\n]+',j);
		idlist=re.findall('[^;]+',elements[17]);
		plist=re.findall('[^;]+',elements[4]);
		#print(idlist);print(plist);
		for k in range(len(idlist)):
			tmp=re.findall('[^-]+',idlist[k]);
			id=elements[1]+'-'+tmp[0]+'-'+tmp[1];
			if i == 0:
				this_inc1=tmp[0]+'_'+tmp[1];
				this_skp1=tmp[0]+'_'+tmp[1];
				id2p[id]="{:.2E}".format(Decimal(float(plist[k])))
				#print('this_inc1');print(this_inc1);print('\n');
			if id in id2p:
				if this_inc1 in Inc:
					Inc[this_inc1].append(type_list[i]+'_'+id+'_'+id2p[id]);
				else:
					Inc[this_inc1]=[type_list[i]+'_'+id+'_'+id2p[id]];
				# if this_inc2 in Inc:
					# Inc[this_inc2].append(type_list[i]+'_'+id+'_'+id2p[id]);
				# else:
					# Inc[this_inc2]=[type_list[i]+'_'+id+'_'+id2p[id]];
				if this_skp1 in Skp:
					Skp[this_skp1].append(type_list[i]+'_'+id+'_'+id2p[id]);
				else:
					Skp[this_skp1]=[type_list[i]+'_'+id+'_'+id2p[id]];
				# if this_skp2 in Skp:
					# Skp[this_skp2].append(type_list[i]+'_'+id+'_'+id2p[id]);
				# else:
					# Skp[this_skp2]=[type_list[i]+'_'+id+'_'+id2p[id]];
	ifile.close();
	
#read the junctions for each ASM
ifile=open(sys.argv[2]);
ilines=ifile.readlines();
ofile=open(sys.argv[3],'w');
ofile.write('ASMID\tMatch\tSE\tA3SS\tA5SS\tMXE\tSEFDR5\tA3SSFDR5\tA5SSFDR5\tMXEFDR5\tPmin\n')
IncEvent = []; SkpEvent=[]; thisASM = '';
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
			p = [0,0,0,0]; pmax = 0;
			for j in thisEvent:
				thisoutput += j+',';
				elements2 = re.findall('[^_]+',j);
				for k in range(len(type_list)):
					if elements2[0] == type_list[k]:
						pmax = max(pmax, float(elements2[2]))
						#print(float(elements2[2]));
						#print(min(pmin, float(elements2[2])));
						#print('===');
						if float(elements2[2]) <= FDRCut[k]:
							p[k] += 1;
			poutput = '';
			for j in p:
				poutput += '\t'+str(j);							
			if len(thisEvent) > 0:
				thisoutput = thisoutput[:-1];
			else:
				thisoutput = 'NA';
			ofile.write(thisoutput+eventoutput+poutput+'\t'+str(pmax)+'\n');
		#write information for the current ASM
		thisASM = elements[0];	iscoord = 1;
		IncEvent = []; SkpEvent=[];
		continue;
	#read coordinates
	if iscoord == 1:
		iscoord = 0; juc1 = []; juc2 = [];
		for j in elements:
			jucs = re.findall('[^,]+',j);
			juc1.append(str(int(jucs[0])));
			juc2.append(jucs[1]);
	#read junctions
	else:
		for j in range(len(elements)-2):
			thisjuc=juc2[int(elements[j+1])]+'_'+juc1[int(elements[j+2])];
			if thisjuc in Inc:
				IncEvent+=Inc[thisjuc];
			if thisjuc in Skp:
				SkpEvent+=Skp[thisjuc];
