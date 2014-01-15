x = read.table(file = "./R_in", header =F,sep="\t")
x[,3]=x[,2]-x[,1];
mean_len = mean(x[,3]);
write.table(mean_len, file ="R_out",quote = F, sep ="\t", row.names = F, col.names = F);
q("no");