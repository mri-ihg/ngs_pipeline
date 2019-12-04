library("RMySQL")

###########################################
#
# This script calculates stretches of homozygosity for a given sample and inserts them into a database.
# It takes the required variant information from the same database.
#
# Parameters are: username, password, database name and host and the name of the sample
#
###########################################

############################# get parameters & connect db #################################
myargs   <-commandArgs(trailingOnly = TRUE)
user     <-myargs[1]
password <-myargs[2]
dbname   <-myargs[3]
host     <-myargs[4]
name     <-myargs[5]
coredb   <-myargs[6]

con <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)

############################# sliding window #################################

# http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r
# sliding window

slideSum<-function(x,windowsize=5,slide=1){
 idx1<-seq(1,length(x),by=slide);
 idx1+windowsize->idx2;
 idx2[idx2>(length(x)+1)]<-length(x)+1;
 c(0,cumsum(x))->cx;
 return((cx[idx2]-cx[idx1]));
}
########################### intodb #################################

intodb<-function(idsample,chromdb,startdb,enddb,n)
{
	SQL<-paste("INSERT IGNORE into homozygosity (idhomozygosity,idsample,chrom,start,end,count) VALUES (NULL,",idsample,",'",chromdb,"',",startdb,",",enddb,",",n,");",sep="")
	dbGetQuery(con, SQL)
}
############################# main searches #################################

# search for idsample
SQL<-paste(
"SELECT idsample
FROM ",coredb,".sample
WHERE name ='",name,"'",
sep='')
idsample<-dbGetQuery(con, SQL)

# search number of samples finished
SQL<-"
SELECT
count(distinct idsample)
FROM snvsample"
nsamples<-dbGetQuery(con, SQL)

#Replace old joins with tmp tables with less items
#INNER JOIN snvsample x ON v.idsnv=x.idsnv
#INNER JOIN ",coredb,".sample s ON x.idsample = s.idsample
# search for alleles
SQL <- paste(
"SELECT v.chrom,v.start,x.alleles,1-freq/",nsamples,
"FROM snv v
INNER JOIN (SELECT * from snvsample where idsample =", idsample, " ) x ON v.idsnv=x.idsnv
INNER JOIN (SELECT * from ",coredb,".sample where idsample=", idsample, ") s ON x.idsample = s.idsample
WHERE s.idsample = ", idsample ,"
AND x.mapqual >= 50
AND x.snvqual >= if(x.caller='samtools',30,3)
AND x.coverage >= 20
AND v.freq >= 100
AND v.class = 'snp'
AND x.alleles<=2
ORDER BY v.chrom,v.start
")
alldata<-dbGetQuery(con, SQL)

############################ delete from database before insert ##########################
SQL<-paste("DELETE FROM homozygosity WHERE idsample=",idsample,sep="")
dbGetQuery(con, SQL)


############################# score #################################

# Wang, Ott et al
# observed  autozygous        not_autozygous
# AA        (1-e)*pa+e*pa^2    pa^2
# AB         2*e*pa*pb         2*pa*pb
# BB        (1-e)*pb+e*pb^2    pb^2


# pa allele frequency of allele a
# pb allele frequency of allele b
# e refers to the rate of genotyping errors and mutations
# input column1 genotype AA=0 AB=1 BB=2
# input column2 allele frequency of a
# Pemperton window size = 60

e<-0.001
ic<-1    # increment
w<- 41   # window always odd

chroms<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
'chr20','chr21','chr22')

for (chrom in chroms) {

sl<-alldata[,1] == chrom
mydata<-alldata[sl,]

# if genotype AA=0
sl<-mydata[,3] == 0
mydata[sl,5]<-log10(((1-e)*mydata[sl,4]+e*mydata[sl,4]^2)/(mydata[sl,4]^2))

# if genotype AB=1
sl<-mydata[,3] == 1
mydata[sl,5]<-log10((2*e*mydata[sl,4]*(1-mydata[sl,4]))/(2*mydata[sl,4]*(1-mydata[sl,4])))

# if genotype BB=2
sl<-mydata[,3] == 2
mydata[sl,5]<-log10(((1-e)*(1-mydata[sl,4])+e*(1-mydata[sl,4])^2)/((1-mydata[sl,4])^2))


#for (i in c(1:(l-w+1))) {
#mydata[(i+(w-1)/2),6]<-sum(mydata[i:(i+w-1),5])
#}

mydata[,6]<-slideSum(mydata[,5],w,ic)
mydata[,7]<-slideSum(mydata[,5],w,ic)





l<-length(mydata[,1])
start<-1
n<-0
chromdb<-""
startdb<-""
enddb<-""
for (i in c(1:l)) {

	# homozygous region start
	if (mydata[i,5]>=0.1 && mydata[i,7]>=3 && start) {
		cat(mydata[i,1])
		chromdb<-mydata[i,1]
		cat(" ")
		cat(mydata[i,2])
		startdb<-mydata[i,2]
		cat(" ")
		start<-0
		n<-0
	}
	n<-n+1
	# homozygous region end
	if (mydata[i,5]<=0.1 && mydata[i,7]<=0 && start==0) {
		cat(mydata[i,1])
		cat(" ")
		cat(mydata[i,2])
		enddb<-mydata[i,2]
		cat(" n=")
		cat(n)
		cat(" i=")
		cat(i)
		cat("\n")
		start<-1
		intodb(idsample,chromdb,startdb,enddb,n)
	}
	# chromosome end
	if (i==l && start==0) {
		cat(mydata[i,1])
		cat(" ")
		cat(mydata[i,2])
		enddb<-mydata[i,2]
		cat(" n=")
		cat(n)
		cat(" i=")
		cat(i)
		cat("\n")
		start<-1
		intodb(idsample,chromdb,startdb,enddb,n)
	}
}
} #for

############################################################
test<-function()
{
plot(density(mydata[,6],bw=1,na.rm=T))
hist(mydata[,7],breaks=100)


sl<-mydata[,1]=='chr6'
test<-mydata[sl,]
plot(test[,2],test[,7],type='l')
plot(density(test[,6],bw=1,na.rm=T))
hist(test[,7],breaks=100)
test[300:600,]
test[1:100,]
mydata[32954:33342,]
}
