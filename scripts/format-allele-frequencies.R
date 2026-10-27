# Check scaffold lengths (to set up a different run, start on line 17)
fai = read.table("/home/centos/USS/erik/lemurs/source_files/Varecia_rubra.fa.fai",header=F)
head(fai)
barplot(fai$V2)

# Combine the allele frequencies with the chromosomal positions
freqs = read.table("/home/centos/USS/erik/lemurs/source_files/Varu_Andranobe_2019.frq",header=T)
pos = read.table("/home/centos/USS/erik/lemurs/source_files/AllVarecia_SUPERautoz.pos",header=F)
names(pos) = c("CHR","POS")
dat = cbind(freqs,pos$POS)
head(dat)
names(dat)[7] = "POS"
dat = dat[,c("CHR","POS","MAF")]
dat = dat[dat$MAF > 0,]

write.table(dat,"/home/centos/USS/erik/lemurs/source_files/Varu_Andranobe_2019_withPosition.frq",quote=F,row.names=F)
dat = read.table("/home/centos/USS/erik/lemurs/source_files/Varu_Andranobe_2019_withPosition.frq",header=T)

# Extract only the freqs from the chromosome of interest, in this case SUPER_12
chr = na.omit(dat[dat$CHR=="SUPER_11",])

# Finally, given an input sample size (census size), convert these frequencies into HWE genotypes by random assignment
K = 150
genos = matrix("0/0",nrow=nrow(chr),ncol=K)
for(i in 1:nrow(chr)){
  freq = chr$MAF[i];
  het_count = round((freq*(1-freq)*2)*K)
  hom_count = round((freq^2)*K)
  
  targets = sort(sample((1:K),(het_count+hom_count)))
  hom_targets = sort(sample(targets,hom_count))
  het_targets = setdiff(targets,hom_targets)
  genos[i,het_targets] = "0/1"
  genos[i,hom_targets] = "1/1"
}
genos = as.data.frame(genos)
names(genos) = paste0("i",1:K)

# Add prefix columns to make this vcf
prefix = data.frame(matrix(".",nrow=nrow(chr),ncol=9))
names(prefix) = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
prefix$`#CHROM`=1
prefix$POS = chr$POS
prefix$REF = "A"
prefix$ALT = "T"
prefix$QUAL = "1000"
prefix$FORMAT = "GT"

# And join
output = "/home/centos/USS/erik/lemurs/Vava_Ranomafana_simulations/input_files/sim-genome-setup_SUPER11_N150.vcf"
vcf = cbind(prefix,genos)
write.table("##fileformat=VCFv4.2",output,row.names=F,col.names=F,quote=F)
write.table(vcf,output,row.names=F,quote=F,append=T,sep="\t")
