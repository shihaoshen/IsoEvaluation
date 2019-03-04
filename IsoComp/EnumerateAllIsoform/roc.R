library(AUC)


plot_roc=function(x,true,color){
	naindex = (!is.na(x)) & (!is.na(true))
    #print(length(x));print(length(true));print(sum(naindex));
	if (sum(naindex) > 5){
        x=x[naindex];
        true=true[naindex];
        #x=x[order(true, decreasing = TRUE)]
        #true=sort(true,decreasing = TRUE);
        true=true[order(x)];
        #print(true[1:10]);
        tmp <- true == 1
        tmp[!is.na(true)] <- cumsum(tmp[!is.na(true)])
        tpr = tmp/sum(true==1,na.rm=T)

        tmp <- true == 0
        tmp[!is.na(true)] <- cumsum(tmp[!is.na(true)])
        fpr = tmp/sum(true==0,na.rm=T)
        
        lines(y=tpr,x=fpr,type='l',xlab='False Positive Rate',ylab='True Positive Rate',col=color);
        
        true_na = true[!is.na(true)]
        pred = seq(length(true_na),1,-1)/length(true_na)
        thisroc = roc(pred,as.factor(true_na))
        print(auc(thisroc))
        
        #calculate tpr at 1% and 5% fpr
        tpr1 = thisroc$tpr[which.min(abs(thisroc$fpr - 0.01))];
        tpr5 = thisroc$tpr[which.min(abs(thisroc$fpr - 0.05))];
        tpr10 = thisroc$tpr[which.min(abs(thisroc$fpr - 0.1))];
        
        return(trunc(c(auc(thisroc), tpr1, tpr5, tpr10)*10^2)/10^2)
    }
}

addlabel = function(label, auc, tpr1, tpr5, tpr10){
    newlabel=NULL;
    for (i in 1:length(label)){
        newlabel=c(newlabel, paste(label[i],auc[i],tpr1[i],tpr5[i],tpr10[i],sep=' | '))
    }
    return(newlabel);
}

#roc for iso
data=read.table('simu_res_10isocount_v3_glmmTrue.txt',header=T,sep='\t')
true = rep(NA,length(data[,7]))
true[abs(data[,7]) > 0.1]=1
true[abs(data[,7]) < 0.01]=0
column=c(6,8,9,10,11);
color=c('black','red','blue','green','brown');
title=' | AUC | TPR at 1% FPR | 5% | 10%'
label=c('rMATS-ISO','rMATS-turbo','majiq','leafcutter','jum');
pdf('rMATS_ISO_Flux_roc.pdf');
plot(x=c(0,1),y=c(0,1),type='n',xlab='False positive rate',ylab='True positive rate',main='Gold Standard Filtered by Total Count >=10');

auc=NULL;tpr1=NULL;tpr5=NULL;tpr10=NULL;
for (i in 1:5){
    if (i %in% c(1:2,4:5)){d = data[,column[i]];
        }else{d = 1 - data[,column[i]];}
	tmp = plot_roc(d,true,color[i]);
    auc=c(auc, tmp[1]); tpr1=c(tpr1,tmp[2]); tpr5=c(tpr5,tmp[3]); tpr10=c(tpr10,tmp[4]);
}
newlabel = addlabel(label, auc, tpr1, tpr5, tpr10);
legend(x=0.3,y=0.3,legend=c(title,newlabel),col=c('white',color),lty=1);

#SE event only
auc=NULL;tpr1=NULL;tpr5=NULL;tpr10=NULL;
plot(x=c(0,1),y=c(0,1),type='n',xlab='False positive rate',ylab='True positive rate',main='SE only');
for (i in 1:5){
    if (i %in% c(1:2,4:5)){d = data[,column[i]];
        }else{d = 1 - data[,column[i]];}
	tmp = plot_roc(d[data[,12]=='SE'],true[data[,12]=='SE'],color[i]);
    auc=c(auc, tmp[1]); tpr1=c(tpr1,tmp[2]); tpr5=c(tpr5,tmp[3]); tpr10=c(tpr10,tmp[4]);
}
newlabel = addlabel(label, auc, tpr1, tpr5, tpr10);
legend(x=0.3,y=0.3,legend=c(title,newlabel),col=c('white',color),lty=1);


#ASS event only
auc=NULL;tpr1=NULL;tpr5=NULL;tpr10=NULL;
plot(x=c(0,1),y=c(0,1),type='n',xlab='False positive rate',ylab='True positive rate',main='ASS only');
for (i in 1:5){
    if (i %in% c(1:2,4:5)){d = data[,column[i]];
        }else{d = 1 - data[,column[i]];}
	tmp=plot_roc(d[data[,12]=='ASS'],true[data[,12]=='ASS'],color[i]);
    auc=c(auc, tmp[1]); tpr1=c(tpr1,tmp[2]); tpr5=c(tpr5,tmp[3]); tpr10=c(tpr10,tmp[4]);
}
newlabel = addlabel(label, auc, tpr1, tpr5, tpr10);
legend(x=0.3,y=0.3,legend=c(title,newlabel),col=c('white',color),lty=1);


#SE + ASS event only
auc=NULL;tpr1=NULL;tpr5=NULL;tpr10=NULL;
plot(x=c(0,1),y=c(0,1),type='n',xlab='False positive rate',ylab='True positive rate',main='SE and ASS');
for (i in 1:5){
    if (i %in% c(1:2,4:5)){d = data[,column[i]];
        }else{d = 1 - data[,column[i]];}
	tmp=plot_roc(d[data[,12]=='ASS' | data[,12]=='SE'],true[data[,12]=='ASS' | data[,12]=='SE'],color[i]);
    auc=c(auc, tmp[1]); tpr1=c(tpr1,tmp[2]); tpr5=c(tpr5,tmp[3]); tpr10=c(tpr10,tmp[4]);
}
newlabel = addlabel(label, auc, tpr1, tpr5, tpr10);
legend(x=0.3,y=0.3,legend=c(title,newlabel),col=c('white',color),lty=1);

dev.off()
