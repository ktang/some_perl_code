x = read.table(file = "./R_in", header =T,sep="\t")
cor_val = cor(x);
#print(cor_val);
x[,3]=x[,1]-x[,2];
mean_diff = mean(x[,3]);
median_diff = median(x[,3]);
#print(mean_diff);
write.table(cor_val,file ="out.txt",quote = F,sep= "\t",row.names=F,col.names =F);
write.table(mean_diff,file ="out.txt",append = T, quote = F,sep= "\t",row.names=F,col.names =F);
write.table(median_diff,file ="out.txt",append =T,quote = F,sep= "\t",row.names=F,col.names =F);
q("no");