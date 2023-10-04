#!/usr/bin/R
#########################################################
#  identify CNS 
# ====================================================
# by Zhai Hang, 2023
#########################################################


# colliner gene calculate #########################################
We used MScsan in tbtools software for the computation of 
co-lineage genes in maize and other grasses.
# genome：
maize B73 (AGPv4), Zea mexican (Yang et al., 2017) , 
sorghum (Sorghum bicolor) (Tao et al., 2021), foxtail millet (Setaria italica), 
two different coix lacryma-jobi L. and perennial grass (Guo et al., 2020). 
###################################################################


# Getting non-coded regions using R ###############################
suppressPackageStartupMessages({library(tidyverse);
  library(data.table);
  library(ggpubr);
  library(ggsci);
  library(doParallel);
  library(foreach);
  library(openxlsx)})

options(scipen=200)


### 1.Counting colinear gene sets in maize and other plants ####
# 1.1 Obtaining data on co-linear genes
BJcoix = "~/Synteny/maize_BJcoix/maizev4.faBJCoxi_v1.fa.geneLinks.tab.xls";
Lachesis = "~/Synteny/maize_Lachesis/maizev4.faLachesis.fa.geneLinks.tab.xls";
Setaria_italica = "~/Synteny/maize_Setaria_italica/maizev4.faSetaria_italica_v2.0.fa.geneLinks.tab.xls";
Sorghum_bicolor = "~/analysis/Synteny/maize_Sorghum_bicolor/maizev4.faSorghum_bicolor_v3.fna.geneLinks.tab.xls";
Mexicana = "~/analysis/Synteny/maize_Mexicana/maizev4.faMexicana.fa.geneLinks.tab.xls";


getDat <- function(x,y){
    if (x in c("BJcoix","Lachesis")){
        tmp <- read.table(x,header = F,stringsAsFactors = F,sep = '\t')
        colnames(tmp) = c("maize",x)
        which(nchar(tmp$maize) == 11)->index
        tmp[nchar(tmp$maize) == y,c(1,2)] = tmp[nchar(tmp$maize) == y,c(2,1)]
    }


    if (x in c("Setaria_italica","bicolor")){
      tmp <- read.table(x,header = F,stringsAsFactors = F,sep = '\t')
      colnames(tmp) = c("maize",x)  
    }

    
    if (x in c("Mexicana")){
      tmp <- read.table(x,header = F,stringsAsFactors = F,sep = '\t')
      tmp$V1 = gsub(".+-(.+)","\\1",tmp$V1)
      tmp$V2 = gsub(".+-(.+)","\\1",tmp$V2)
      colnames(tmp) = c("maize",x)  
    }
    return(tmp)
}
getDat(x = "BJcoix",y = 11) -> BJcoix
getDat(x = "Lachesis",y = 12) -> Lachesis
getDat(x = "Setaria_italica") -> Setaria_italica
getDat(x = "bicolor") -> bicolor
getDat(x = "Mexicana") -> Mexicana

raw_syntenic.df <- merge(merge(merge(merge(BJcoix,Lachesis,by = "maize",all = T),
                                    Setaria_italica,by = "maize",all = T),
                               Sorghum_bicolor,by = "maize",all = T),
                        Mexicana,by = "maize",all = T)

# maize
read.table("~/Genomes/maize_v4/maizev4.gtf",header=F,stringsAsFactors=F,sep='\t')->maize.gtf
maize.gtf = maize.gtf[maize.gtf$V3 == "transcript",]
maize.gtf$gid = gsub("gene_id ([^;]+); transcript_id ([^;]+);.+","\\1",maize.gtf$V9)
maize.gtf$tid = gsub("gene_id ([^;]+); transcript_id ([^;]+);.+","\\2",maize.gtf$V9)

# BJCoxi
read.table("~/Genomes/BJCoxi_v1/BJCoix_v1.gff3",header=F,stringsAsFactors=F,sep='\t')->BJcoix.gff3
BJcoix.gff3 = BJcoix.gff3[BJcoix.gff3$V3 == "mRNA",]
BJcoix.gff3$gid = gsub("ID=([^;]+)_.+;.+","\\1",BJcoix.gff3$V9)
BJcoix.gff3$tid = gsub("ID=([^;]+);.+","\\1",BJcoix.gff3$V9)

# Lachesis
read.table("~/Genomes/Lachesis/Lachesis.gff3",header=F,stringsAsFactors=F,sep='\t')->Lachesis.gff3
Lachesis.gff3 = Lachesis.gff3[Lachesis.gff3$V3 == "mRNA",]
Lachesis.gff3$gid = gsub("ID=([^;]+)\\..+;.+","\\1",Lachesis.gff3$V9)
Lachesis.gff3$tid = gsub("ID=([^;]+);.+","\\1",Lachesis.gff3$V9)

# Sorghum_bicolor
read.table("~/Genomes/Sorghum_bicolor_v3/Sorghum_bicolor_v3.gff",header=F,stringsAsFactors=F,sep='\t')->Sorghum_bicolor.gff
Sorghum_bicolor.gff = Sorghum_bicolor.gff[Sorghum_bicolor.gff$V3 == "mRNA",]
Sorghum_bicolor.gff$gid = gsub(".+Parent=(.+);","\\1",Sorghum_bicolor.gff$V9)
Sorghum_bicolor.gff$tid = gsub("ID=([^;]+);.+","\\1",Sorghum_bicolor.gff$V9)

# Mexicana
read.table("~/Genomes/Mexicana/Mexicana.gff3",header=F,stringsAsFactors=F,sep='\t')->Mexicana.gff3
Mexicana.gff3 = Mexicana.gff3[Mexicana.gff3$V3 == "mRNA",]
Mexicana.gff3$gid = gsub(".+Parent=([^;].+);.+","\\1",Mexicana.gff3$V9)
Mexicana.gff3$tid = gsub(".+ID=(.+)","\\1",Mexicana.gff3$V9)

# Setaria_italica
read.table("~/Genomes/Setaria_italica_v2.0/Setaria_italica_v2.0.gff",header=F,quote = "$",stringsAsFactors=F,sep='\t')->Setaria_italica_v2.0.gff
Setaria_italica_v2.0.gff = Setaria_italica_v2.0.gff[Setaria_italica_v2.0.gff$V3 == "mRNA",]
Setaria_italica_v2.0.gff$gid = gsub(".+gene=([^;]+);.+","\\1",Setaria_italica_v2.0.gff$V9)
Setaria_italica_v2.0.gff$tid = gsub(".+transcript_id=(.+)","\\1",Setaria_italica_v2.0.gff$V9)

syntenic.df = raw_syntenic.df
syntenic.df$maize = maize.gtf[match(syntenic.df$maize,maize.gtf$tid),"gid"]
syntenic.df$BJcoix = BJcoix.gff3[match(syntenic.df$BJcoix,BJcoix.gff3$tid),"gid"]
syntenic.df$Lachesis = Lachesis.gff3[match(syntenic.df$Lachesis,Lachesis.gff3$tid),"gid"]
syntenic.df$Sorghum_bicolor = Sorghum_bicolor.gff[match(syntenic.df$Sorghum_bicolor,Sorghum_bicolor.gff$tid),"gid"]
syntenic.df$Mexicana = Mexicana.gff3[match(syntenic.df$Mexicana,Mexicana.gff3$tid),"gid"]
syntenic.df$Setaria_italica = Setaria_italica_v2.0.gff[match(syntenic.df$Setaria_italica,Setaria_italica_v2.0.gff$tid),"gid"]

# 1.2 Number of colinear genes 
syntenic.df$number = apply(syntenic.df,1,function(x){
  length(which(!is.na(x)))
})
syntenic.df$number = syntenic.df$number - 1

# 2.Gene upstream and downstream 10kb bed files
maize = "~/Genomes/maize_v4/maizev4.gtf";
BJcoix = "~/Genomes/BJCoxi_v1/BJCoix_v1.gff3";
Lachesis = "~/Genomes/Lachesis/Lachesis.gff3";
Sorghum_bicolor = "~/Genomes/Sorghum_bicolor_v3/Sorghum_bicolor_v3.gff";
Mexicana = "~/Genomes/Mexicana/Mexicana.gff3";
Setaria_italica = "~/Genomes/Setaria_italica_v2.0/Setaria_italica_v2.0.gff";

range10kb <- function(x){
    read.table(x,header=F,stringsAsFactors=F,sep='\t') -> tmp
    tmp = tmp[tmp$V3 == "gene",]

# maize
    if (x == "maize"){
        tmp$gid = gsub("gene_id ([^;]+);.+","\\1",tmp$V9)
    }
# BJcoix
    if (x == "BJcoix"){
        tmp$gid = gsub("ID=([^;]+);Name.+","\\1",tmp$V9)
    }
# Lachesis
    if (x == "Lachesis"){
        tmp$gid = gsub("ID=([^;]+);","\\1",tmp$V9)
    }
# Sorghum_bicolor
    if (x == "Sorghum_bicolor"){
        tmp$gid = gsub("ID=([^;]+);","\\1",tmp$V9)
    }
# Mexicana    
    if (x == "Mexicana"){
        tmp$gid = gsub(".+ID=([^;]+).+","\\1",tmp$V9)  
    }
# Setaria_italica
    if (x == "Setaria_italica"){
        tmp$gid = gsub("ID=gene-([^;]+);.+","\\1",tmp$V9)
    }


tmp.bed = data.frame(chr = tmp$V1,
                        start = tmp$V4,
                        end = tmp$V5,
                        id = tmp$gid,
                        score = ".",
                        strand = tmp$V7,
                        stringsAsFactors = F)
tmp_10k.bed = tmp.bed
tmp_10k.bed$new_start = tmp_10k.bed$start - 10000
tmp_10k.bed$new_end = tmp_10k.bed$end + 10000
tmp_10k.bed$new_start[tmp_10k.bed$new_start <=0] = 1 
tmp_10k.bed = tmp_10k.bed[,c(1,7,8,4,5,6,2,3)]

return(tmp.10k.bed)
}

range10kb(x = "maize") -> maize_10k.bed
range10kb(x = "BJcoix") -> BJcoix_10k.bed
range10kb(x = "Lachesis") -> Lachesis_10k.bed
range10kb(x = "Sorghum_bicolor") -> Sorghum_bicolor_10k.bed
range10kb(x = "Mexicana") -> Mexicana_10k.bed
range10kb(x = "Setaria_italica") -> Setaria_italica_10k.bed


# 3.Co-linear genes upstream and downstream 10kb bed files
# 3.1 maize
maize_10k.syntenic.bed = maize_10k.bed[maize_10k.bed$id %in% syntenic.df$maize,]
write.table(maize_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > maize_10k.syntenic.bed")

write.table(maize.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > maize.bed")

# 3.2 BJcoix 
BJcoix_10k.syntenic.bed = BJcoix_10k.bed[BJcoix_10k.bed$id %in% syntenic.df$BJcoix,]
write.table(BJcoix_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > BJcoix_10k.syntenic.bed")

write.table(BJcoix.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > BJcoix.bed")

# 3.3 Lachesis
Lachesis_10k.syntenic.bed = Lachesis_10k.bed[Lachesis_10k.bed$id %in% syntenic.df$Lachesis,]
write.table(Lachesis_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Lachesis_10k.syntenic.bed")

write.table(Lachesis.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Lachesis.bed")

# 3.4 Sorghum_bicolor
Sorghum_bicolor_10k.syntenic.bed = Sorghum_bicolor_10k.bed[Sorghum_bicolor_10k.bed$id %in% syntenic.df$Sorghum_bicolor,]
write.table(Sorghum_bicolor_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Sorghum_bicolor_10k.syntenic.bed")

write.table(Sorghum_bicolor.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Sorghum_bicolor.bed")

# 3.5 Mexicana
Mexicana_10k.syntenic.bed = Mexicana_10k.bed[Mexicana_10k.bed$id %in% syntenic.df$Mexicana,]
write.table(Mexicana_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Mexicana_10k.syntenic.bed")

write.table(Mexicana.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Mexicana.bed")

# 3.6 Setaria_italica
Setaria_italica_10k.syntenic.bed = Setaria_italica_10k.bed[Setaria_italica_10k.bed$id %in% syntenic.df$Setaria_italica,]
write.table(Setaria_italica_10k.syntenic.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Setaria_italica_10k.syntenic.bed")

write.table(Setaria_italica.bed,"tmp.bed",col.names = F,row.names = F,quote = F,sep = '\t')
system("bedtools sort -i tmp.bed > Setaria_italica.bed")


#### 4.Non-coding regions between covariate genes:
# Check if there is a co-linear gene within a 10-kilobase range upstream and downstream of a co-linear gene.
# If such a co-linear gene exists, extract the non-coding sequence of that gene located between another co-linear gene.
# If the gene in question is not a non-co-linear gene, do not make any changes.
### function1: Detect other co-linear genes within the upstream and downstream range.
cal_syntenic.fun1 <- function(input.bed.r,input.bed,ref.bed.r,ref.bed,mat){
  
  system(paste("bedtools intersect -a ",input.bed.r," -b ",ref.bed.r," -wo >tmp.out",sep = ""))
  tmp.out <- read.table("tmp.out",header = F,stringsAsFactors = F,sep = '\t')
  tmp.out = tmp.out[tmp.out$V4 != tmp.out$V12,] 
  
  ## # Retrieve other co-linear genes within the upstream and downstream regions that overlap with the co-linear gene.
  tmp.out = tmp.out[tmp.out$V12 %in% syntenic.df[,mat],]
  
  ddply(tmp.out,.(V4),.parallel = T,function(x){
    filename=tempfile(tmpdir=".") 
    tmp1.bed = paste(filename,"tmp1.bed",sep = "")
    tmp2.bed = paste(filename,"tmp2.bed",sep = "")
    
    write.table(input.bed[input.bed$id == unique(x$V4),1:6],tmp1.bed,row.names = F,col.names = F,quote = F,sep = '\t')#本底
    write.table(ref.bed[ref.bed$id %in% x$V12,1:6],tmp2.bed,row.names = F,col.names = F,quote = F,sep = '\t')#其它非线性基因
    
    system(paste("bedtools subtract -a ",tmp1.bed," -b ",tmp2.bed,">",filename,sep = "")) #取非overlap区域
    readLines(filename)->tmp
    
    if (length(tmp) > 0) {
      read.table(filename,header = F,stringsAsFactors = F,sep = '\t')->tmp
      write.table(tmp[,1:6],tmp2.bed,row.names = F,col.names = F,quote = F,sep = '\t')#非overlap区域
      write.table(ref.bed[ref.bed$id == unique(x$V4),1:6],tmp1.bed,row.names = F,col.names = F,quote = F,sep = '\t')#本地基因位点
      
      system(paste("bedtools closest -a ",tmp1.bed," -b ",tmp2.bed,">",filename,sep = "")) #取非overlap区域离本地最近的区域
      read.table(filename,header = F,stringsAsFactors = F,sep = '\t')->tmp
      
      out.bed = tmp[1,7:12]
      out.bed
    }else{
      out.bed = data.frame(V7=NA,V8=NA,V9=NA,V10=unique(x$V4),V11=NA,V12=NA,stringsAsFactors = F) #如果没有非overlap区域，记录下这个位点，后面则舍弃这个位点
      out.bed 
    }
    unlink(filename)
    unlink(tmp1.bed)
    unlink(tmp2.bed)
    out.bed
  })[,-1] -> diff.bed
  
  index = match(diff.bed$V10,input.bed$id)
  tmp = input.bed
  tmp[index,1:6] = diff.bed[,1:6]
  tmp = tmp[!is.na(tmp[,1]),]
  tmp
}
### function2: Organize and process the results.
cal_syntenic.fun2 <- function(input.bed,syntenic.bed){
  tmp.bed = get(input.bed)
  syntenic.bed = get(syntenic.bed)
  ddply(tmp.bed,.(id),.parallel = T,function(x){
    if (unique(x$type == 3)) {
      data.frame(chr = unique(x$chr),
                 new_start = x[1,2],
                 new_end = x[2,3],
                 id = unique(x$id),
                 score = ".",
                 strand = unique(x$strand),
                 start = x[3,2],
                 end = x[3,3],
                 stringsAsFactors = F)->tmp
    }else{
      if (x$type[1] == "up.out") {
        data.frame(chr = unique(x$chr),
                   new_start = x[1,2],
                   new_end = syntenic.bed[syntenic.bed$id == unique(x$id),"new_end"],
                   id = unique(x$id),
                   score = ".",
                   strand = unique(x$strand),
                   start = x[2,2],
                   end = x[2,3],
                   stringsAsFactors = F)->tmp
        
      }else{
        data.frame(chr = unique(x$chr),
                   new_start = syntenic.bed[syntenic.bed$id == unique(x$id),"new_start"],
                   new_end = x[1,3],
                   id = unique(x$id),
                   score = ".",
                   strand = unique(x$strand),
                   start = x[2,2],
                   end = x[2,3],
                   stringsAsFactors = F)->tmp
        
      }
      
    }
    return(tmp)
  })}

### 4.1 maize ####
cal_syntenic.fun1(input.bed.r = "maize_10k.syntenic.bed",
                 input.bed = maize_10k.syntenic.bed,
                 ref.bed.r = "maize.bed",
                 ref.bed = maize.bed,
                 mat = "maize")->tmp.bed

tmp1.bed = maize_10k.syntenic.bed[!maize_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed


tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "maize_10k.syntenic.bed")->tmp2.bed
maize_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)

maize_10k.syntenic_final.bed = rbind(maize_10k.syntenic.bed[!maize_10k.syntenic.bed$id %in% maize_10k.syntenic_final.bed$id,],maize_10k.syntenic_final.bed)

table(maize_10k.syntenic_final.bed$chr) 
maize_10k.syntenic_final.bed = maize_10k.syntenic_final.bed[maize_10k.syntenic_final.bed$chr %in% 1:10,]
write.table(maize_10k.syntenic_final.bed,"maize_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

### 4.2 BJcoix ####
rm(tmp.bed,tmp1.bed,tmp2.bed)
cal_syntenic.fun1(input.bed.r = "BJcoix_10k.syntenic.bed",
                  input.bed = BJcoix_10k.syntenic.bed,
                  ref.bed.r = "BJcoix.bed",
                  ref.bed = BJcoix.bed,
                  mat = "BJcoix")->tmp.bed

tmp1.bed = BJcoix_10k.syntenic.bed[!BJcoix_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed

tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "BJcoix_10k.syntenic.bed")->tmp2.bed
BJcoix_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)
BJcoix_10k.syntenic_final.bed = rbind(BJcoix_10k.syntenic.bed[!BJcoix_10k.syntenic.bed$id %in% BJcoix_10k.syntenic_final.bed$id,],BJcoix_10k.syntenic_final.bed)

table(BJcoix_10k.syntenic_final.bed$chr) #10chr,过滤
BJcoix_10k.syntenic_final.bed = BJcoix_10k.syntenic_final.bed[BJcoix_10k.syntenic_final.bed$chr %in% paste("Chr",1:10,sep = ""),]
write.table(BJcoix_10k.syntenic_final.bed,"BJcoix_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

### 4.3 Lachesis ####
rm(tmp.bed,tmp1.bed,tmp2.bed)

Lachesis_10k.syntenic.bed[Lachesis_10k.syntenic.bed$id == "EVM0010711",]
cal_syntenic.fun1(input.bed.r = "Lachesis_10k.syntenic.bed",
                 input.bed = Lachesis_10k.syntenic.bed,
                 ref.bed.r = "Lachesis.bed",
                 ref.bed = Lachesis.bed,
                 mat = "Lachesis")->tmp.bed

tmp1.bed = Lachesis_10k.syntenic.bed[!Lachesis_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed

tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "Lachesis_10k.syntenic.bed")->tmp2.bed
Lachesis_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)
Lachesis_10k.syntenic_final.bed = rbind(Lachesis_10k.syntenic.bed[!Lachesis_10k.syntenic.bed$id %in% Lachesis_10k.syntenic_final.bed$id,],Lachesis_10k.syntenic_final.bed)

table(Lachesis_10k.syntenic_final.bed$chr)
Lachesis_10k.syntenic_final.bed = Lachesis_10k.syntenic_final.bed[Lachesis_10k.syntenic_final.bed$chr %in% paste("Lachesis_group0",0:9,sep = ""),]
write.table(Lachesis_10k.syntenic_final.bed,"Lachesis_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

### 4.4 Sorghum_bicolor ####
rm(tmp.bed,tmp1.bed,tmp2.bed)

cal_syntenic.fun1(input.bed.r = "Sorghum_bicolor_10k.syntenic.bed",
                  input.bed = Sorghum_bicolor_10k.syntenic.bed,
                  ref.bed.r = "Sorghum_bicolor.bed",
                  ref.bed = Sorghum_bicolor.bed,
                  mat = "Sorghum_bicolor")->tmp.bed

tmp1.bed = Sorghum_bicolor_10k.syntenic.bed[!Sorghum_bicolor_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed

tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "Sorghum_bicolor_10k.syntenic.bed")->tmp2.bed
Sorghum_bicolor_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)
Sorghum_bicolor_10k.syntenic_final.bed = rbind(Sorghum_bicolor_10k.syntenic.bed[!Sorghum_bicolor_10k.syntenic.bed$id %in% Sorghum_bicolor_10k.syntenic_final.bed$id,],Sorghum_bicolor_10k.syntenic_final.bed)

table(Sorghum_bicolor_10k.syntenic_final.bed$chr) #10chr,无需过滤
write.table(Sorghum_bicolor_10k.syntenic_final.bed,"Sorghum_bicolor_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

### 4.5 Mexicana ####
rm(tmp.bed,tmp1.bed,tmp2.bed)

cal_syntenic.fun1(input.bed.r = "Mexicana_10k.syntenic.bed",
                  input.bed = Mexicana_10k.syntenic.bed,
                  ref.bed.r = "Mexicana.bed",
                  ref.bed = Mexicana.bed,
                  mat = "Mexicana")->tmp.bed

tmp1.bed = Mexicana_10k.syntenic.bed[!Mexicana_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed

tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "Mexicana_10k.syntenic.bed")->tmp2.bed
Mexicana_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)
Mexicana_10k.syntenic_final.bed = rbind(Mexicana_10k.syntenic.bed[!Mexicana_10k.syntenic.bed$id %in% Mexicana_10k.syntenic_final.bed$id,],Mexicana_10k.syntenic_final.bed)

table(Mexicana_10k.syntenic_final.bed$chr) #10chr,无需过滤
write.table(Mexicana_10k.syntenic_final.bed,"Mexicana_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

### 4.6 Setaria_italica ####
rm(tmp.bed,tmp1.bed,tmp2.bed)

cal_syntenic.fun1(input.bed.r = "Setaria_italica_10k.syntenic.bed",
                  input.bed = Setaria_italica_10k.syntenic.bed,
                  ref.bed.r = "Setaria_italica.bed",
                  ref.bed = Setaria_italica.bed,
                  mat = "Setaria_italica")->tmp.bed

tmp1.bed = Setaria_italica_10k.syntenic.bed[!Setaria_italica_10k.syntenic.bed$id %in% tmp.bed$id,]

ddply(tmp.bed,.(id),.parallel = T,function(x){
  x$num = nrow(x)
  
  if (x$num == 3) {
    x$log_type = all(c(x[1,3] == x[3,2],x[2,2] == x[3,3]))
    
  }else{
    if (x$type[1] == "up.out") {
      x$log_type = x[1,3] == x[2,2]
    }else{
      x$log_type = x[1,2] == x[2,3]
    }
  }
  x
})->tmp.bed

tmp.bed = tmp.bed[tmp.bed$log_type == TRUE,]
tmp.bed = tmp.bed[!is.na(tmp.bed$start),]

cal_syntenic.fun2(input.bed = "tmp.bed",syntenic.bed = "Setaria_italica_10k.syntenic.bed")->tmp2.bed
Setaria_italica_10k.syntenic_final.bed = rbind(tmp1.bed,tmp2.bed)
Setaria_italica_10k.syntenic_final.bed = rbind(Setaria_italica_10k.syntenic.bed[!Setaria_italica_10k.syntenic.bed$id %in% Setaria_italica_10k.syntenic_final.bed$id,],Setaria_italica_10k.syntenic_final.bed)

table(Setaria_italica_10k.syntenic_final.bed$chr) #10chr,无需过滤
write.table(Setaria_italica_10k.syntenic_final.bed,"Setaria_italica_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

#### 5. filter ####
## 5.1 Filter length ( < 8bp )
## maize
tmp <- read.table("~/Genomes/maize_v4/maizev4.fa.fai",header = F,stringsAsFactors = F,sep = '\t')[,1:2]
ddply(maize_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->tmp1
  maize_10k.syntenic_final.bed
write.table(tmp1,"maize_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## BJcoix
tmp <- read.table("~/Genomes/BJCoxi_v1/BJCoxi_v1.fa.fai",header = F,stringsAsFactors = F,sep = '\t')[1:10,1:2]
ddply(BJcoix_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->BJcoix_10k.syntenic_final.bed

write.table(BJcoix_10k.syntenic_final.bed,"BJcoix_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Lachesis
tmp <- read.table("~/Genomes/Lachesis/Lachesis.fa.fai",header = F,stringsAsFactors = F,sep = '\t')[1:10,1:2]
ddply(Lachesis_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->Lachesis_10k.syntenic_final.bed
write.table(Lachesis_10k.syntenic_final.bed,"Lachesis_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Sorghum_bicolor
tmp <- read.table("~/Genomes/Sorghum_bicolor_v3/Sorghum_bicolor_v3.fna.fai",header = F,stringsAsFactors = F,sep = '\t')[1:10,1:2]
ddply(Sorghum_bicolor_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->Sorghum_bicolor_10k.syntenic_final.bed

write.table(Sorghum_bicolor_10k.syntenic_final.bed,"Sorghum_bicolor_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Mexicana
tmp <- read.table("~/Genomes/Mexicana/Mexicana.fa.fai",header = F,stringsAsFactors = F,sep = '\t')[1:10,1:2]
tmp$V1 = paste("chr",tmp$V1,sep = "")
ddply(Mexicana_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->Mexicana_10k.syntenic_final.bed

write.table(Mexicana_10k.syntenic_final.bed,"Mexicana_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Setaria_italica_10k.syntenic_final.bed
tmp <- read.table("~/Genomes/Setaria_italica_v2.0/Setaria_italica_v2.0.fa.fai",header = F,stringsAsFactors = F,sep = '\t')
tmp = tmp[grep("NC",tmp$V1),]
ddply(Setaria_italica_10k.syntenic_final.bed,.(chr),.parallel = T,function(x){
  x = x[x$new_start < tmp[tmp$V1 == unique(x$chr),2],]
  x$new_end[x$new_end >= tmp[tmp$V1 == unique(x$chr),2]] = tmp[tmp$V1 == unique(x$chr),2]
  x$len = x$new_end - x$new_start
  x = x[x$len > 8,]
})->Setaria_italica_10k.syntenic_final.bed

write.table(Setaria_italica_10k.syntenic_final.bed,"Setaria_italica_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## 5.2 Elimination of unpaired genes
## maize
ddply(maize_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->maize_10k.syntenic_final.bed

write.table(maize_10k.syntenic_final.bed,"maize_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## BJcoix
ddply(BJcoix_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->BJcoix_10k.syntenic_final.bed

write.table(BJcoix_10k.syntenic_final.bed,"BJcoix_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Lachesis
ddply(Lachesis_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->Lachesis_10k.syntenic_final.bed

write.table(Lachesis_10k.syntenic_final.bed,"Lachesis_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Sorghum_bicolor
ddply(Sorghum_bicolor_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->Sorghum_bicolor_10k.syntenic_final.bed

write.table(Sorghum_bicolor_10k.syntenic_final.bed,"Sorghum_bicolor_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Mexicana
ddply(Mexicana_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->Mexicana_10k.syntenic_final.bed
Mexicana_10k.syntenic_final.bed$chr = gsub("chr","",Mexicana_10k.syntenic_final.bed$chr)
write.table(Mexicana_10k.syntenic_final.bed,"Mexicana_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

## Setaria_italica_10k.syntenic_final.bed
ddply(Setaria_italica_10k.syntenic_final.bed,.(gid),.parallel = T,function(x){
  x$pair = nrow(x)
  x = x[x$pair == 2,]
  x
})->Setaria_italica_10k.syntenic_final.bed

write.table(Setaria_italica_10k.syntenic_final.bed,"Setaria_italica_10k.syntenic_final.bed",col.names = F,row.names = F,quote = F,sep = '\t')

#### 6. Preparing to input data
## 6.1 extract sequence
library(Biostrings)
system("bedtools getfasta -fi ~/Genomes/maize_v4/maizev4.fa -name -bed maize_10k.syntenic_final.bed>maize_10k.syntenic_final.fa")
system("bedtools getfasta -fi ~/Genomes/BJCoxi_v1/BJCoxi_v1.fa -name -bed BJcoix_10k.syntenic_final.bed>BJcoix_10k.syntenic_final.fa")
system("bedtools getfasta -fi ~/Genomes/Lachesis/Lachesis.fa -name -bed Lachesis_10k.syntenic_final.bed>Lachesis_10k.syntenic_final.fa")
system("bedtools getfasta -fi ~/Genomes/Sorghum_bicolor_v3/Sorghum_bicolor_v3.fna -name -bed Sorghum_bicolor_10k.syntenic_final.bed>Sorghum_bicolor_10k.syntenic_final.fa")
system("bedtools getfasta -fi ~/Genomes/Mexicana/Mexicana.fa -name -bed Mexicana_10k.syntenic_final.bed>Mexicana_10k.syntenic_final.fa")
system("bedtools getfasta -fi ~/Genomes/Setaria_italica_v2.0/Setaria_italica_v2.0.fa -name -bed Setaria_italica_10k.syntenic_final.bed>Setaria_italica_10k.syntenic.fa")

## 6.2 get.type function: extract non-coding sequence 
get.type = function(input.fa,input.bed){
  readDNAStringSet(input.fa)->fa
  data.frame(seq = as.character(fa),
             name = names(fa),
             stringsAsFactors = F) ->seq
  seq$gid = gsub("(.+)::.+","\\1",seq$name)
  seq$chr = gsub(".+::(.+):(.+)-(.+)","\\1",seq$name)
  seq$start = gsub(".+::(.+):(.+)-(.+)","\\2",seq$name)
  seq$end = gsub(".+::(.+):(.+)-(.+)","\\3",seq$name)
  seq$strand = input.bed[match(seq$gid,input.bed$id),"strand"]
  
  for (i in colnames(seq)[-c(1,2)]) {
    print(i)
    print(table(nchar(seq[,i])))
  } 
  return(seq)
}

get.type("maize_10k.syntenic_final.fa",input.bed = maize_10k.syntenic.bed)->maize_10k.syntenic_final.seq
get.type("BJcoix_10k.syntenic_final.fa",input.bed = BJcoix_10k.syntenic.bed)->BJcoix_10k.syntenic_final.seq
get.type("Lachesis_10k.syntenic_final.fa",input.bed = Lachesis_10k.syntenic.bed)->Lachesis_10k.syntenic_final.seq
get.type("Sorghum_bicolor_10k.syntenic_final.fa",input.bed = Sorghum_bicolor_10k.syntenic.bed)->Sorghum_bicolor_10k.syntenic_final.seq
get.type("Mexicana_10k.syntenic_final.fa",input.bed = Mexicana_10k.syntenic.bed)->Mexicana_10k.syntenic_final.seq
get.type("Setaria_italica_10k.syntenic.fa",input.bed = Setaria_italica_10k.syntenic.bed)->Setaria_italica_10k.syntenic_final.seq

### 6.3 re-filter
syntenic_clean.df = syntenic.df

## maize
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$maize),"maize"] %in% maize_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$maize)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]

## BJcoix
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$BJcoix),"BJcoix"] %in% BJcoix_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$BJcoix)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]

## Lachesis
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$Lachesis),"Lachesis"] %in% Lachesis_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$Lachesis)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]

## Sorghum_bicolor
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$Sorghum_bicolor),"Sorghum_bicolor"] %in% Sorghum_bicolor_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$Sorghum_bicolor)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]

## Mexicana
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$Mexicana),"Mexicana"] %in% Mexicana_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$Mexicana)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]

## Setaria_italica
index1 = syntenic_clean.df[!is.na(syntenic_clean.df$Setaria_italica),"Setaria_italica"] %in% Setaria_italica_10k.syntenic_final.seq$gid
index2 = rownames(syntenic_clean.df)[!is.na(syntenic_clean.df$Setaria_italica)][!index1]

syntenic_clean.df = syntenic_clean.df[!row.names(syntenic_clean.df) %in% index2,]


## count
for (i in colnames(syntenic_clean.df)) {
  print(i)
  print(length(unique(syntenic_clean.df[,i])))
  print(table(is.na(syntenic_clean.df[,i])))
} 

for (i in colnames(syntenic.df)) {
  print(i)
  print(length(unique(syntenic.df[,i])))
  print(table(is.na(syntenic.df[,i])))
} 

nrow(syntenic_clean.df)
nrow(syntenic.df)
table(syntenic_clean.df$number)

### 6.4 Input file for STAG-CNS
data.frame(index = 1:nrow(syntenic_clean.df),tmp = "-",stringsAsFactors = F)->index
ddply(index,.(index),.parallel = T,function(x){
  
  maize = maize_10k.syntenic_final.seq
  BJcoix = BJcoix_10k.syntenic_final.seq
  Lachesis = Lachesis_10k.syntenic_final.seq
  Sorghum_bicolor = Sorghum_bicolor_10k.syntenic_final.seq
  Mexicana = Mexicana_10k.syntenic_final.seq
  Setaria_italica = Setaria_italica_10k.syntenic_final.seq
  
  tmp = syntenic_clean.df[x$index,]
  name.index = colnames(tmp)[1:6][which(!is.na(unlist(tmp,)[1:6]))]
  
  info1 = data.frame(tmp = "-")
  info2 = data.frame(tmp = "-")
  
  for (i in 1:length(name.index)) {
    tmp1 = get(name.index[i])
    tmp2 = get(paste(name.index[i],"_10k.syntenic_final.bed",sep = ""))
    id1 = tmp[,name.index[i]]
    
    data.frame(paste(">",id1,
                     unique(tmp1[tmp1$gid == id1,"chr"]),
                     unique(tmp1[tmp1$gid == id1,"strand"]),
                     unique(tmp1[tmp1$gid == id1,"start"]),
                     unique(tmp2[tmp2$id == id1,"start"]),
                     unique(tmp2[tmp2$id == id1,"end"]),
                     unique(tmp1[tmp1$gid == id1,"end"]),
                     sep = "\t"))->paste1
    info1 = cbind(info1,paste1)
    
    data.frame(paste(">",id1,sep = ""),
               tmp1[tmp1$gid == id1,"seq"])->paste2
    info2 = cbind(info2,paste2)
  }
  
  info1 = info1[,-1]
  info2 = info2[,-1]
  info = cbind(info1,info2)

  write.table(info,paste("new_input/",x$index,".",length(name.index),"_input.txt",sep = ""),col.names = F,row.names = F,quote = F,sep = "\n")
  x$index
})->tmp


#### 6.5 Compute in Cluster
#!/bin/sh
#PBS -l nodes=16:ppn=4
#PBS -j oe
#PBS -o run.log

# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
cd $PBS_O_WORKDIR

cat ../analysis/cmat | /opt/parallel/bin/parallel -j8 --sshloginfile $PBS_NODEFILE --workdir $PWD --env PATH --env LD_LIBRARY_PATH  '
cd ../analysis/new_output
mkdir {}_output
cd {}_output

stagCNS -file ../../new_input/{}_input.txt -mem 8 -out {} > {}.out.txt
'

