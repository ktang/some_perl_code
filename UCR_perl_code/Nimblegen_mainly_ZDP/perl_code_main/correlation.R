x = read.table(file = "./R_in", header =T,sep="\t")
cor_val = cor(x);
print(cor_val);
x[,3]=x[,1]-x[,2];
mean_diff = mean(x[,3]);
print(mean_diff);
q("no");