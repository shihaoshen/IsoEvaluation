data=read.table('glmer_input.txt',header=T);
data$group = as.factor(data$group)

#try glmm first
data$psi = data$inc / data$tot;
#print(data$psi); print(data$group);
#res = try(wilcox.test(psi ~ group, data = data))
res = try(t.test(psi ~ group, data = data))

#glmm not works, use glm
if (attr(res,'class')!='try-error'){
	p = res$p.value;
	#print(res);
}else{
	p = 1;
}
write.table(p,file='ttest_output.txt');
