library(AUC)


plot_roc=function(x,true,color){
	naindex = (!is.na(x)) & (!is.na(true))
	x=x[naindex];
	true=true[naindex];
	#x=x[order(true, decreasing = TRUE)]
	#true=sort(true,decreasing = TRUE);
	true=true[order(x)];
	print(true[1:10]);
	tmp <- true == 1
	tmp[!is.na(true)] <- cumsum(tmp[!is.na(true)])
	tpr = tmp/sum(true==1,na.rm=T)

	tmp <- true == 0
	tmp[!is.na(true)] <- cumsum(tmp[!is.na(true)])
	fpr = tmp/sum(true==0,na.rm=T)
	
	lines(y=tpr,x=fpr,type='l',xlab='False Positive Rate',ylab='True Positive Rate',col=color);
	
	true_na = true[!is.na(true)]
	pred = seq(length(true_na),1,-1)/length(true_na)
	roc(pred,as.factor(true_na))
	print(auc(roc(pred,as.factor(true_na))))
}
	
#roc for iso
data=read.table('simu_res.txt',header=T,sep='\t')
true = rep(NA,length(data[,7]))
true[abs(data[,7]) > 0.1]=1
true[abs(data[,7]) < 0.01]=0
column=c(6,8,9,10,11);
color=c('black','red','blue','green','brown');
label=c('ISO','turbo','majiq','leafcutter','jum');
pdf('rMATS_ISO_Flux_roc.pdf');
plot(x=c(0,1),y=c(0,1),type='n',xlab='False positive rate',ylab='True positive rate')

for (i in 1:5){
	plot_roc(data[,column[i]],true,color[i]);
}
legend(x=0.6,y=0.3,legend=label,col=color,lty=1);
dev.off()


