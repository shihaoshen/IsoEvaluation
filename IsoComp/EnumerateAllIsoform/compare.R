#this script compares the 
library(AUC)

#roc for iso
data=read.table('ASM2ASMdetails.txt',header=T,sep='\t')

#in SE
p_new = data$iso_p_value;
p_old = data$iso_p_value.1;
type = data$Type
na_index=(!is.na(p_new)) & (!is.na(p_old)) & (type == 'SE' | type == 'ASS') & (!is.na(type))
print('Number of matched iso');
sum(na_index);
print('Number of matched iso with identical p');
sum((p_new[na_index] - p_old[na_index])<10^-5)
print('Number of matched iso with p < 0.05 in both cases');
sum((p_new[na_index] < 0.05) & (p_old[na_index]<0.05))
print('Number of matched iso with p < 0.05 in one cases');
sum((p_new[na_index] < 0.05) & (p_old[na_index]>0.05)) + sum((p_new[na_index] > 0.05) & (p_old[na_index]<0.05))
pdf('rMATS_ISO_P_with_isoforms_enumeration_vs_major_isoform_search_SEonly.pdf')
plot(y=p_new[na_index],x=p_old[na_index],ylab='rMATS-ISO p with all isoforms outputed',xlab='rMATS-ISO p with major isoform search before statistical model',main='rMATS-ISO p value in SE & ASS')
plot(y=rank(p_new[na_index]),x=rank(p_old[na_index]),ylab='rMATS-ISO p rank with all isoforms',xlab='rMATS-ISO p rank with major isoform search before statistical model',main='rMATS-ISO p value rank')
dev.off()

#in all modules
p_new = data$iso_p_value;
p_old = data$iso_p_value.1;
na_index=(!is.na(p_new)) & (!is.na(p_old))
print('Number of matched iso');
sum(na_index);
print('Number of matched iso with identical p');
sum((p_new[na_index] - p_old[na_index])<10^-5)
print('Number of matched iso with p < 0.05 in both cases');
sum((p_new[na_index] < 0.05) & (p_old[na_index]<0.05))
print('Number of matched iso with p < 0.05 in one cases');
sum((p_new[na_index] < 0.05) & (p_old[na_index]>0.05)) + sum((p_new[na_index] > 0.05) & (p_old[na_index]<0.05))
pdf('rMATS_ISO_P_with_isoforms_enumeration_vs_major_isoform_search_allisoforms.pdf')
all=sum(na_index);
equal=sum((p_new[na_index] - p_old[na_index])<10^-5)
plot(y=p_new[na_index],x=p_old[na_index],ylab='rMATS-ISO p with all isoforms outputed',xlab='rMATS-ISO p with major isoform search before statistical model',main=paste('rMATS-ISO p value in all modules\n',equal,' out of ',all,' modules have identical p values',sep=''))
plot(y=rank(p_new[na_index]),x=rank(p_old[na_index]),ylab='rMATS-ISO p rank with all isoforms',xlab='rMATS-ISO p rank with major isoform search before statistical model',main='rMATS-ISO p value rank')
dev.off();


