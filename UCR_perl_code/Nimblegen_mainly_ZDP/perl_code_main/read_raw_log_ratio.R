setwd();

files = dir();

raw = NULL

##################

raw = read.table (file = files[1],header=F,sep ="\t");
raw = raw[,3];
raw = as.data.frame(raw);
colnames(raw) = files[1];

for (i in 2:3)
{
	temp = read.table (file = files[i],header=F,sep ="\t")
	temp = temp[,3]
	temp = as.data.frame(temp);
	colnames(temp) = files[i];
	raw[,i] = temp[,1];
}

for (i in 4:10)
{
	temp = read.table (file = files[i],header=F,sep ="\t")
	temp = temp[,3];
	temp = as.data.frame(temp);
	colnames(temp) = files[i];
	raw[,i] = temp[,1];
}
############

raw = read.table (file = files[1],header=F,sep ="\t");
raw = raw[,3];
raw = as.data.frame(raw);

for (i in 2:10)
{
	temp = read.table (file = files[i],header=F,sep ="\t")
	temp = temp[,3]
	temp = as.data.frame(temp);

	raw[,i] = temp;
}

colnames(raw) = files;

for (i in 1:9)
{
	for (j in (i+1) : 10)
	{
		print (  c (cor(raw[,i],raw[,j]), colnames(raw)[i], colnames(raw)[j]))
	}
}

for (i in 1:9)
{
	for (j in (i+1) : 10)
	{
		write.table  (  paste (cor(raw[,i],raw[,j]), colnames(raw)[i], colnames(raw)[j] ,sep="\t")  , file ="cor_for_all_45.txt", append = T,quote=F,sep="",
					  row.names = F,col.names=F)
	}
}


for (i in 1:9)
{
	for (j in (i+1) : 10)
	{
		print (  c(mean(raw[,i]-raw[,j]), colnames(raw)[i], colnames(raw)[j]))
	}
}

for (i in 1:9)
{
	for (j in (i+1) : 10)
	{
		print (  paste(mean(raw[,i]-raw[,j]), colnames(raw)[i], colnames(raw)[j]) ,sep ="\t")
	}
}



for (i in 1:9)
{
	for (j in (i+1) : 10)
	{
		write.table  (  paste (mean(raw[,i]-raw[,j]), colnames(raw)[i], colnames(raw)[j] ,sep="\t")  , file ="mean_diff_for_all_45.txt", append = T,quote=F,sep="",
					  row.names = F,col.names=F)
	}
}
