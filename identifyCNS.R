#!/usr/bin/R
#########################################################
#  identify CNS 
# ====================================================
# by Zhai Hang, 2023
#########################################################


# colliner gene calculate #########################################
We used MScsan in tbtools software for the computation of 
co-lineage genes in maize and other grasses.
(maize B73 (AGPv4), Zea mexican (Yang et al., 2017) , 
sorghum (Sorghum bicolor) (Tao et al., 2021), foxtail millet (Setaria italica), 
two different coix lacryma-jobi L. and perennial grass (Guo et al., 2020)). 
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
