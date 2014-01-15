setwd("");


#############
for (i in 1:1)
{
	c = RGs_separate[i];
	log_ratio = log2(c[[1]]$R) - log2(c[[1]]$G);
	rownames(log_ratio) = c[[1]]$genes$ID;
	write.table(log_ratio,file = colnames(log_ratio),quote = F,sep ="\t",
				col.names = F);
}




for (i in 2:10)
{
	c = RGs_separate[i];
	log_ratio = log2(c[[1]]$R) - log2(c[[1]]$G);
	rownames(log_ratio) = c[[1]]$genes$ID;
	write.table(log_ratio,file = colnames(log_ratio),quote = F,sep ="\t",
				col.names = F);
}
###################



for (i in 1:10)
{
	c = RGs_separate[i];
	log_ratio = log2(c[[1]]$R) - log2(c[[1]]$G);
	rownames(log_ratio) = c[[1]]$genes$ID;
	write.table(log_ratio,file = colnames(log_ratio),quote = F,sep ="\t",
				col.names = F);
}