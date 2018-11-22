#Requires genes to be in the first column

library(tidyverse)
library(broom)
library(shiny)
library(GenABEL)
library(readxl)

#Genes names in multiple formats
gm <- read_tsv("Data/Gene_names/gene_names.txt")

#Get 1to1 orthologs from Jax database
m2h <- read_tsv("Data/Mouse2human/m2h.txt")

#Function to read file and process

read.files.process <- function(df){
  cell <- read.table(df,header=T,stringsAsFactors = FALSE)
  colnames(cell) <- c("musName","t","P","fdr")
  cell <- as.tibble(cell)
  cell_1to1 <- dplyr::inner_join(cell,m2h,by="musName")
  cell_1to1 <- dplyr::select(cell_1to1,geneName,t,P)
  colnames(cell_1to1) <- c("Gene","t","pvalue")
  return(cell_1to1)
}

files <- list.files(path = "Data/DE/",pattern=".BPSC.txt",full.names = TRUE)

d <- data_frame(filename = files) %>% mutate(file_contents = purrr::map(filename,read.files.process)) %>%
  unnest() %>% mutate(Cell_type=gsub("\\."," ",gsub(".+\\/\\/","",gsub("\\.vs.all.BPSC.txt","",filename)))) %>% 
  dplyr::filter(!is.na(t)) %>% group_by(Cell_type) %>% mutate(t_norm=rntransform(t))

#View(d %>% group_by(Cell_type) %>% count())

#tbl <- read.csv("~/Desktop/Gene_sets_random/extTADA.txt",stringsAsFactors = FALSE, sep="\t")
tbl <- read_xlsx("~/Desktop/Gene_sets_random/tada.top288.xlsx")
colnames(tbl)[1] <- "Gene"
  
#If name start with ENSG strip possible numerical values after .
if (all(grepl("ENSG00",tbl$Gene))){
  tbl <- mutate(tbl,Gene=gsub("\\..+","",tbl$Gene))
}
    
#Select column with the most match to the input gene list
column_to_merge <- which.max(c(length(which(tbl$Gene%in%gm$gene_name)),
                  length(which(tbl$Gene%in%gm$entrez_id)),
                  length(which(tbl$Gene%in%gm$ensembl_id)),
                  length(which(tbl$Gene%in%gm$hgnc_id)),
                  length(which(tbl$Gene%in%gm$refseq_accession))))
    
#If the first column has the most matches, keep original IDs that were tested in the DE (genes with a 1to1 ortholog)
#else merge the gene_list with the gene name annotations to get the official gene name
    
if(column_to_merge==1){
  tbl <- tbl %>% filter(Gene %in% d$Gene)
      
  } else{
    tbl <- inner_join(tbl,gm,by=c("Gene" = colnames(gm)[column_to_merge]))
    tbl <- tbl[-1]
    tbl <- rename(tbl,Gene=gene_name)
    tbl <- tbl %>% dplyr::filter(Gene %in% d$Gene)
    is.numeric <- sapply(tbl, is.numeric)
    tbl <- cbind(tbl["Gene"],tbl[is.numeric])
}
  
#minimum <- 0
#maximum <- 0.1
#filtered_df <- tbl[na.omit(tbl$qvalue_DD >= minimum & tbl$qvalue_DD<=maximum),]

#d <- d %>% mutate(gene_set=ifelse(Gene%in%filtered_df$Gene,1,0),all_genes=ifelse(Gene%in%tbl$Gene,1,0))
#d <- d %>% group_by(Cell_type) %>% mutate(t_norm=rntransform(t))
#regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ all_genes + gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="gene_set") %>% arrange(P)
#regression_log <- d %>% group_by(Cell_type) %>% do(tidy(glm(gene_set ~ all_genes + t_norm,family=binomial(link='logit'),data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="t_norm") %>% arrange(P)
#wilcox_results <- d %>% group_by(Cell_type) %>% do(tidy(wilcox.test(.$t_norm[.$gene_set==1],.$t_norm[.$gene_set==0],alternative="greater"))) %>% arrange(p.value)

#regression_t <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ all_genes + gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="gene_set") %>% arrange(P)
#regression_log_t <- d %>% group_by(Cell_type) %>% do(tidy(glm(gene_set ~ all_genes + t,family=binomial(link='logit'),data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="t") %>% arrange(P)
#wilcox_results_t <- d %>% group_by(Cell_type) %>% do(tidy(wilcox.test(.$t[.$gene_set==1],.$t[.$gene_set==0],alternative="greater"))) %>% arrange(p.value)

filtered_df <- tbl
d <- d %>% mutate(gene_set=ifelse(Gene%in%filtered_df$Gene,1,0))
d <- d %>% group_by(Cell_type) %>% mutate(t_norm=rntransform(t))
regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="gene_set") %>% arrange(P)

#Plot example
d <- d %>% ungroup %>% mutate(Cell_type=factor(Cell_type,levels=regression$Cell_type))
vlines <- d %>% group_by(Cell_type,gene_set) %>% summarise(median=median(t_norm),mean=mean(t_norm)) %>% mutate(col=ifelse(gene_set==0,"#67a9cf","#ef8a62"))

p <- d  %>% ggplot(aes(t_norm,color=gene_set)) + geom_histogram(bins=100) + facet_wrap(~Cell_type)
p <- p + theme_bw() + xlab("Standard normalised t statistic") + theme(strip.text = element_text(size=7))
ggsave(p,filename = "QC/t_norm/Plots/hist_per_cell_type_all_genes_eli_288_genes.pdf",width=12)

p1 <- ggplot(d,aes(t_norm))
p1 <- p1 + geom_histogram(data=subset(d,gene_set == 1),fill = "#ef8a62", alpha = 0.5,aes(y=..count../sum(..count..)))
p1 <- p1 +  geom_histogram(data=subset(d,gene_set == 0),fill = "#67a9cf", alpha = 0.5,aes(y=..count../sum(..count..))) 
p1 <- p1 + theme_bw() + theme(strip.text = element_text(size=7))
p1 <- p1 + geom_vline(data=vlines,aes(xintercept=mean,color=col),alpha=1)
p1 <- p1 + facet_wrap(~Cell_type) + scale_colour_identity() + xlab("Standard normalised t statistic") + ylab("Frequency")
ggsave(p1,filename = "QC/t_norm/Plots/hist_per_cell_type_eli_288_genes.pdf",width=12)

p1 <- ggplot(d,aes(t_norm))
p1 <- p1 + geom_histogram(data=subset(d,gene_set == 1),fill = "#ef8a62", alpha = 0.5,aes(y=..density..))
p1 <- p1 +  geom_histogram(data=subset(d,gene_set == 0),fill = "#67a9cf", alpha = 0.5,aes(y=..density..))
p1 <- p1 + theme_bw() + theme(strip.text = element_text(size=7))
p1 <- p1 + geom_vline(data=vlines,aes(xintercept=mean,color=col),alpha=1)
p1 <- p1 + facet_wrap(~Cell_type) + scale_colour_identity() + xlab("Standard normalised t statistic") + ylab("Frequency")
ggsave(p1,filename = "QC/t_norm/Plots/hist_per_cell_type_eli_288_genes2.pdf",width=12)

p1 <- ggplot(d,aes(t_norm))
p1 <- p1 + geom_freqpoly(data=subset(d,gene_set == 1),col = "#ef8a62", alpha = 0.5,aes(y=..density..))
p1 <- p1 + geom_freqpoly(data=subset(d,gene_set == 0),col = "#67a9cf", alpha = 0.5,aes(y=..density..))
p1 <- p1 + theme_bw() + theme(strip.text = element_text(size=7))
p1 <- p1 + geom_vline(data=vlines,aes(xintercept=mean,color=col),alpha=1)
p1 <- p1 + facet_wrap(~Cell_type) + scale_colour_identity() + xlab("Standard normalised t statistic") + ylab("Frequency")
ggsave(p1,filename = "QC/t_norm/Plots/hist_per_cell_type_eli_288_genes3.pdf",width=12)

p1 <- ggplot(d,aes(t_norm))
p1 <- p1 + geom_density(data=subset(d,gene_set == 1),col = "#ef8a62", alpha = 0.2,fill="#ef8a62")
p1 <- p1 + geom_density(data=subset(d,gene_set == 0),col = "#67a9cf", alpha = 0.2,fill="#67a9cf")
p1 <- p1 + theme_bw() + theme(strip.text = element_text(size=7))
p1 <- p1 + geom_vline(data=vlines,aes(xintercept=mean,color=col),alpha=1)
p1 <- p1 + facet_wrap(~Cell_type) + scale_colour_identity() + xlab("Standard normalised t statistic") + ylab("Frequency")
ggsave(p1,filename = "QC/t_norm/Plots/hist_per_cell_type_eli_288_genes4.pdf",width=12)



#Test random gene sets Function
    
random_genes <- function(size,n_perm,n_cells){
  random_gene_sets <- matrix(ncol=n_perm,nrow=n_cells)
  for(i in 1:n_perm){
    if(i%%50==0){
      cat(i)
      cat("\n")
    }
    random <- sample(unique(tbl$Gene),size,replace = F)
    d <- d %>% mutate(gene_set=ifelse(Gene%in%random,1,0),all_genes=ifelse(Gene%in%tbl$Gene,1,0))
    regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ all_genes + gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="gene_set")
    random_gene_sets[,i] <- regression$P
  }
  random_gene_sets <- as.data.frame(random_gene_sets)
  rownames(random_gene_sets) <- regression$Cell_type
  write.table(random_gene_sets,paste0("QC/t_norm/random_gene_sets_",n_perm,"_permutations_",n_cells,"_celltypes_",size,"_genes_tnorm.txt"),col.names = F,row.names=T,sep="\t",quote=F)
  return(random_gene_sets)
}    
 
random_50 <- random_genes(50,1000,24)
random_100 <- random_genes(100,1000,24)
random_200 <- random_genes(200,1000,24)
random_300 <- random_genes(300,1000,24)
random_400 <- random_genes(400,1000,24)
random_500 <- random_genes(500,1000,24)
random_1000 <- random_genes(1000,1000,24)
random_2000 <- random_genes(2000,1000,24)
random_3000 <- random_genes(3000,1000,24)
random_5000 <- random_genes(5000,1000,24)
    
#Function to gather all values in table

gather2table <- function(df){
   df_tbl <- df %>% mutate(Cell_type=rownames(df)) %>% gather(key=Cell_type,value=P)
   colnames(df_tbl) <- c("Cell_type","Permutation","P")
   return(df_tbl)
}

random_50_2df <- gather2table(random_50)
random_100_2df <- gather2table(random_100)
random_200_2df <- gather2table(random_200)
random_300_2df <- gather2table(random_300)
random_400_2df <- gather2table(random_400)
random_500_2df <- gather2table(random_500)
random_1000_2df <- gather2table(random_1000)
random_2000_2df <- gather2table(random_2000)
random_3000_2df <- gather2table(random_3000)
random_5000_2df <- gather2table(random_5000)

combined <- rbind(random_50_2df,random_100_2df,random_200_2df,random_300_2df,
                  random_400_2df,random_500_2df,random_1000_2df,
                  random_2000_2df,random_3000_2df,random_5000_2df)

#Function to make qqplot for each cell type

qq <- function(df,title){
df <- df %>% group_by(Cell_type) %>% 
  mutate(random_log10_sort=sort(-log10(seq(1/((nrow(.)/24)),1,1/((nrow(.)/24))))),
         P_log10_sort=sort(-log10(P)))

p1 <- ggplot(df,aes(random_log10_sort,P_log10_sort)) 
p1 <- p1 + geom_point() + facet_wrap(~Cell_type) + geom_abline()
p1 <- p1  + theme_bw() + xlab("Expected (-log10(P))") + ylab("Observed (-log10(P))")
p1 <- p1 + ggtitle(title) + theme(strip.text.x = element_text(size = 5))
p1
}

pdf("QC/t_norm/Plots/Expected_vs_observed.pdf")
qq(random_50_2df,"50 Random genes")
qq(random_100_2df,"100 Random genes")
qq(random_200_2df,"200 Random genes")
qq(random_300_2df,"300 Random genes")
qq(random_400_2df,"400 Random genes")
qq(random_500_2df,"500 Random genes")
qq(random_1000_2df,"1000 Random genes")
qq(random_2000_2df,"2000 Random genes")
qq(random_3000_2df,"3000 Random genes")
qq(random_5000_2df,"5000 Random genes")
qq(combined,"combined Random genes")
dev.off()

### WITH T

random_gene_sets_100 <- read.table("QC/random_gene_sets_100_1K_permutations.txt",sep="\t",header=F,stringsAsFactors = FALSE)
random_gene_sets_200 <- read.table("QC/random_gene_sets_200_1K_permutations.txt",sep="\t",header=F,stringsAsFactors = FALSE)
random_gene_sets_500 <- read.table("QC/random_gene_sets_500_1K_permutations.txt",sep="\t",header=F,stringsAsFactors = FALSE)
random_gene_sets_1000 <- read.table("QC/random_gene_sets_1000_1K_permutations.txt",sep="\t",header=F,stringsAsFactors = FALSE)
random_gene_sets_2000 <- read.table("QC/random_gene_sets_2000_1K_permutations.txt",sep="\t",header=F,stringsAsFactors = FALSE)

gather2table <- function(df){
  df_tbl <- df %>% rename(Cell_type=V1) %>% gather(key=Cell_type,value=P)
  colnames(df_tbl) <- c("Cell_type","Permutation","P")
  return(df_tbl)
}

random_100_2df <- gather2table(random_gene_sets_100)
random_200_2df <- gather2table(random_gene_sets_200)
random_500_2df <- gather2table(random_gene_sets_500)
random_1000_2df <- gather2table(random_gene_sets_1000)
random_2000_2df <- gather2table(random_gene_sets_2000)

combined <- rbind(random_100_2df,random_200_2df,
                  random_500_2df,random_1000_2df,
                  random_2000_2df)

pdf("QC/Plots/Expected_vs_observed.pdf")
qq(random_100_2df,"100 Random genes")
qq(random_200_2df,"200 Random genes")
qq(random_500_2df,"500 Random genes")
qq(random_1000_2df,"1000 Random genes")
qq(random_2000_2df,"2000 Random genes")
qq(combined,"combined Random genes")
dev.off()

## Plot regression results




##
#Test potential covariates
##
cov <- read.table("Data/Covariates/covariates_for_brainrich.txt",header=T,stringsAsFactors = FALSE)
cov <- dplyr::select(cov,-Gene) %>% rename(Gene=Gene_hs) %>% as_tibble()
d <- left_join(d,cov,by="Gene")
    
#regression_gene_length <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ gene_length_hs,data=.))) %>% filter(term=="gene_length_hs")
regression_exon_length <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ exon_length_gencode_v25_hs,data=.))) %>% filter(term=="exon_length_gencode_v25_hs")
regression_sum_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_sum_all,data=.))) %>% filter(term=="Gene_sum_all")
regression_mean_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_mean_all,data=.))) %>% filter(term=="Gene_mean_all")
regression_mean0_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_mean0_all,data=.))) %>% filter(term=="Gene_mean0_all")
    
regression_Gene_var_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_var_all,data=.))) %>% filter(term=="Gene_var_all")
regression_Gene_cv_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_cv_all,data=.))) %>% filter(term=="Gene_cv_all")
    
# Mean + mean0
regression_mean_all_mean0_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_mean_all + Gene_mean0_all,data=.))) %>% filter(term=="Gene_mean0_all")
    
#Mean + mean0 + cv
regression_mean_all_mean0_all_cv_all <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_mean_all + Gene_mean0_all+Gene_cv_all ,data=.))) %>% filter(term=="Gene_cv_all")
    
#Mean + mean0 + exon_length
regression_mean_all_mean0_all_length <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ Gene_mean_all + Gene_mean0_all+exon_length_gencode_v25_hs ,data=.))) %>% filter(term=="exon_length_gencode_v25_hs")
regression2 <- d %>% group_by(Cell_type) %>% do(tidy(lm(t~ all_genes + gene_set + Gene_mean_all + Gene_mean0_all +exon_length_gencode_v25_hs,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="gene_set") %>% arrange(P)