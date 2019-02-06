data=read.table(file='CCLE_EMonly.txt')

cutoff = 5;

for (i in levels(data[,10])){
	EM_stat = summary(data[data[,10]==i,14]);
	if (length(EM_stat)==2){
		if ((EM_stat[1] >= cutoff) & (EM_stat[2] >= cutoff)){
			print(i); print(EM_stat);
		}
	}
}

