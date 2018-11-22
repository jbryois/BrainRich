#Requires genes to be in the first column

# To do
# Add human/Mouse button

library(tidyverse)
library(broom)
library(shiny)
library(readxl)
library(GenABEL)
library(Hmisc)

#options to increase max size of the file to be loaded
options(shiny.maxRequestSize=90*1024^2) 

#Genes names in multiple formats from gene Matrix
gm <- read_tsv("Data/Gene_names/gene_names.txt")

#Get 1to1 orthologs from Jax database
m2h <- read_tsv("Data/Mouse2human/m2h.txt")

#Function to read differential expression files and only keep genes with a 1to1 ortholog

read.files.process <- function(df){
  cell <- read.table(df,header=T,stringsAsFactors = FALSE)
  colnames(cell) <- c("musName","t","P","fdr")
  cell <- as.tibble(cell)
  cell_1to1 <- dplyr::inner_join(cell,m2h,by="musName")
  cell_1to1 <- dplyr::select(cell_1to1,geneName,t,P)
  colnames(cell_1to1) <- c("Gene","t","pvalue")
  #cell_1to1 <- dplyr::select(cell_1to1,geneName,musName,entrez,t,P)
  #colnames(cell_1to1) <- c("Gene","musName","entrez","t","pvalue")
  return(cell_1to1)
}

#Files containing the results of differential expression for each cell type
files <- list.files(path = "Data/DE/",pattern=".BPSC.txt",full.names = TRUE)

d <- data_frame(filename = files) %>% mutate(file_contents = purrr::map(filename,read.files.process)) %>%
  unnest() %>% mutate(Cell_type=gsub("\\."," ",gsub(".+\\/\\/","",gsub("\\.vs.all.BPSC.txt","",filename)))) %>% 
  dplyr::filter(!is.na(t)) %>% mutate(Cell_type=capitalize(Cell_type)) %>% group_by(Cell_type) %>% mutate(t_norm=rntransform(t)) 

# Shiny APP
shinyServer(function(input, output,session) {
  
  #Read the gene-set and change the name of the first column to "Gene"
  gene_list <- reactive({
    inFile <- input$file1
    req(inFile)
    if(!input$sep%in%c('xlsx','xls')){
      tbl <- read.csv(inFile$datapath, header=input$header, sep=input$sep,comment = "#",stringsAsFactors = FALSE)
    } 
    if(input$sep=='xlsx'){
      tbl <- as.data.frame(read_xlsx(inFile$datapath,col_names=input$header))
    }
    if(input$sep=='xls'){
      tbl <- as.data.frame(read_xls(inFile$datapath,col_names=input$header))
    }

    colnames(tbl)[1] <- "Gene"
    return(tbl)
  }) 
  
  #Match first column to gene name ID for multiple different gene name types
  gene_list_id_matched <- reactive({
    
    tbl <- gene_list()
    
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
      tbl <- tbl %>% dplyr::filter(Gene %in% d$Gene)
      is.numeric <- sapply(tbl, is.numeric)
      tbl <- cbind(tbl["Gene"],tbl[is.numeric])
    } else{
      tbl <- inner_join(tbl,gm,by=c("Gene" = colnames(gm)[column_to_merge]))
      tbl <- tbl[-1]
      tbl <- rename(tbl,Gene=gene_name)
      tbl <- tbl %>% dplyr::filter(Gene %in% d$Gene)
      is.numeric <- sapply(tbl, is.numeric)
      tbl <- cbind(tbl["Gene"],tbl[is.numeric])
    }
    return(tbl)
  })
  
  #Sets an UI for selecting a numeric column
  output$choose_column <- renderUI({
    req(ncol(gene_list_id_matched())>1)
    is.numeric <- sapply(gene_list_id_matched(), is.numeric)
    colnames <- colnames(gene_list_id_matched())[is.numeric]
    colnames <- colnames[order(colnames,decreasing = T)]
    req(length(colnames)>=1)
    selectInput("choose_col", "Feature selection:", colnames)
  })
  
  #Sets an UI for selecting a threshold based on the column selected
  output$top <- renderUI({
    req(ncol(gene_list_id_matched())>1)
    req(input$choose_col)
    col_of_interest <- gene_list_id_matched()[,input$choose_col]
    sliderInput("top_n", "Filter:",min=min(col_of_interest,na.rm = T),max=max(col_of_interest,na.rm = T)*1.05,value = c(min(col_of_interest,na.rm = T),quantile(col_of_interest,0.05,na.rm = T)))
    })
  
  #Filter the data frame based on the threshold selected for the column selected
  filter_df <- eventReactive(input$top_n, {
    req(input$top_n)
    minimum <- min(input$top_n)
    maximum <- max(input$top_n)
    filtered_df <- gene_list_id_matched()[na.omit(gene_list_id_matched()[,input$choose_col] >= minimum & gene_list_id_matched()[,input$choose_col]<=maximum),]
    return(filtered_df)
  })
  
  #Plot histogram of the column selected
  output$hist <- renderPlot({
    req(input$choose_col)
    req(ncol(gene_list_id_matched())>1)
    req(input$top_n)
    ggplot(gene_list_id_matched(),aes_string(input$choose_col)) + geom_histogram(bins=100) + theme_bw() + geom_vline(xintercept = input$top_n,color="red")
  })
  
  #Write number of genes matched
  output$Number_genes <- renderTable({
    inFile <- input$file1
    req(inFile)

    #If gene sets as only one column, assume that they are the selected genes
    if (ncol(gene_list_id_matched())==1){
      values <- data.frame(Category=c("Input","Matched","Selected", "Not selected"),
                           Number_genes=c(nrow(gene_list()),nrow(gene_list_id_matched()),nrow(gene_list_id_matched()),"0"))
    } else{
      gene_list_df <- filter_df()
      values <- data.frame(Category=c("Input","Matched","Selected", "Not selected"),
                           Number_genes=c(nrow(gene_list()),nrow(gene_list_id_matched()),nrow(gene_list_df),nrow(gene_list_id_matched())-nrow(gene_list_df)))
    }
    return(values)
  })
  
  reg <- reactive({
    inFile <- input$file1
    req(inFile)
    
    #If gene sets as only one column, assume that they are the selected genes
    if (ncol(gene_list_id_matched())==1){
      req(gene_list_id_matched())
      gene_list_df <- gene_list_id_matched()
      d <- d %>% mutate(gene_set=ifelse(Gene%in%gene_list_df$Gene,1,0))
      regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% dplyr::filter(term=="gene_set") %>% arrange(desc(statistic))
    } 
    
    #For gene list that are small but with multiple columns, the regressions are not corrected for all genes tested.
    if (ncol(gene_list_id_matched())>1 & nrow(gene_list_id_matched())<4000){
      req(filter_df())
      gene_list_df <- filter_df()
      d <- d %>% mutate(gene_set=ifelse(Gene%in%gene_list_df$Gene,1,0))
      regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% dplyr::filter(term=="gene_set") %>% arrange(desc(statistic))
    } 
    
    #For gene lists that are large, assumes that all genes are tested and controls for the genes that were tested
    if (ncol(gene_list_id_matched())>1 & nrow(gene_list_id_matched())>=4000){
      req(filter_df())
      gene_list_df <- filter_df()
      d <- d %>% mutate(gene_set=ifelse(Gene %in% gene_list_df$Gene,1,0),all_genes=ifelse(Gene %in% gene_list_id_matched()$Gene,1,0))
      regression <- d %>% group_by(Cell_type) %>% do(tidy(lm(t_norm~ all_genes + gene_set,data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% dplyr::filter(term=="gene_set") %>% arrange(desc(statistic))
    }
    #regression_log <- d %>% group_by(Cell_type) %>% do(tidy(glm(gene_set ~ all_genes + t,family=binomial(link='logit'),data=.))) %>% mutate(P=ifelse(estimate<0,1-(p.value/2),p.value/2)) %>% filter(term=="t") %>% arrange(P)
    return(regression)
  })
  
  #Plot Cell type enrichment pvalue plot
  output$Enrichment_p <- renderPlot({
    req(reg())
    p3 <- ggplot(reg(),aes(reorder(Cell_type,statistic),-log10(P),fill=reorder(Cell_type,statistic))) + geom_bar(stat="identity") + coord_flip() +theme_bw()+ theme(legend.position = "none",text = element_text(size=18)) 
    p3 <- p3 + ylab(expression('-log'[10]*'(pvalue)')) + xlab("") + geom_hline(yintercept = -log10(0.05/nrow(reg())))
    p3
  })
  
  #Plot regression table (Not used)
  #output$regression_table <- renderDataTable({
  #  out_tbl <- reg() %>% dplyr::select(Cell_type,estimate,std.error,statistic,P)
  #})
  
  #Save BrainRich cell type enrchment pvalue
  output$save_plot <- downloadHandler(
    filename = function() {
      if(ncol(gene_list_id_matched())>1){
        paste(input$file1,".",input$choose_col,".",min(input$top_n),"-",max(input$top_n),".BrainRich", ".pdf", sep="")
      } else{
        paste(input$file1,"BrainRich", "pdf", sep=".")
      }
    },
    content = function(file) {
      p3 <- ggplot(reg(),aes(reorder(Cell_type,-log10(P)),-log10(P),fill=reorder(Cell_type,-log10(P)))) + geom_bar(stat="identity") + coord_flip() +theme_bw()+ theme(legend.position = "none") 
      p3 <- p3 + ylab(expression('-log'[10]*'(pvalue)')) + xlab("") + geom_hline(yintercept = -log10(0.05/nrow(reg())))
      ggsave(file, plot = p3, device = "pdf")
    }
  )
  
  #Save BrainRich cell type histograms
  output$save_plot2 <- downloadHandler(
    filename = function() {
      if(ncol(gene_list_id_matched())>1){
        paste(input$file1,".",input$choose_col,".",min(input$top_n),"-",max(input$top_n),".BrainRich.density", ".pdf", sep="")
      } else{
        paste(input$file1,"BrainRich.density", "pdf", sep=".")
      }
    },
    content = function(file) {
      req(reg())
      if (ncol(gene_list_id_matched())==1){
        req(gene_list_id_matched())
        gene_list_df <- gene_list_id_matched()
      } else{
        req(filter_df())
        gene_list_df <- filter_df()
      }  
      d <- d %>% mutate(gene_set=ifelse(Gene %in% gene_list_df$Gene,1,0))
      d <- d %>% ungroup %>% mutate(Cell_type=factor(Cell_type,levels=reg()$Cell_type))
      vlines <- d %>% group_by(Cell_type,gene_set) %>% summarise(mean=mean(t_norm)) %>% mutate(col=ifelse(gene_set==0,"#67a9cf","#ef8a62"))
      p1 <- ggplot(d,aes(t_norm))
      #p1 <- p1 + geom_histogram(data=subset(d,gene_set == 1),fill = "#ef8a62", alpha = 0.5,aes(y=..count../sum(..count..)))
      #p1 <- p1 + geom_histogram(data=subset(d,gene_set == 0),fill = "#67a9cf", alpha = 0.5,aes(y=..count../sum(..count..))) 
      #p1 <- p1 + geom_density(data=subset(d,gene_set == 1),col = "#ef8a62", alpha = 0.2,fill="#ef8a62")
      #p1 <- p1 + geom_density(data=subset(d,gene_set == 0),col = "#67a9cf", alpha = 0.2,fill="#67a9cf")
      p1 <- p1 + geom_histogram(data=subset(d,gene_set == 1),fill = "#ef8a62", alpha = 0.5,aes(y=..density..),bins=50)
      p1 <- p1 + geom_histogram(data=subset(d,gene_set == 0),fill = "#67a9cf", alpha = 0.5,aes(y=..density..),bins=50) 
      
      p1 <- p1 + theme_bw() + theme(strip.text = element_text(size=7),text = element_text(size=14))
      p1 <- p1 + geom_vline(data=vlines,aes(xintercept=mean,color=col),alpha=1)
      p1 <- p1 + facet_wrap(~Cell_type) + scale_colour_identity() + xlab("Standard normalised t statistic") + ylab("Density")
      ggsave(file, plot = p1, device = "pdf",width=12)
    }
  )
  
  #Save BrainRich table
  output$save_table <- downloadHandler(
    filename = function() {
      if(ncol(gene_list_id_matched())>1){
        paste(input$file1,".",input$choose_col,".",min(input$top_n),"-",max(input$top_n),".BrainRich", ".txt", sep="")
      } else{
        paste(input$file1,"BrainRich", "txt", sep=".")
      }
    },
    content = function(file) {
      out <- reg() %>% dplyr::select(Cell_type,estimate,std.error,statistic,P)
      write_tsv(out, file)
    }
  )

############    
### EWCE  ##
############
  
  ewce_results <- eventReactive(input$ewce_run, {
    load("../BrainRich/Data/Proportion/celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0.rda") 
    proportion <- as.data.frame(celltype_data[[1]]$cell_dists)
    proportion$musName <- rownames(proportion)
    
    #Keep mouse genes with a 1to1 orthoglog and change mouse Name to Human Name
    proportion <- merge(proportion,m2h,by="musName")
    rownames(proportion) <- proportion$geneName
    proportion <- dplyr::select(proportion,astrocytes_ependymal:`Vascular Leptomeningeal Cells`)
    if (ncol(gene_list_id_matched())==1){
      req(gene_list_id_matched())
      gene_list_df <- gene_list_id_matched()
    } else{
      req(filter_df())
      gene_list_df <- filter_df()
    }
    proportion_gene_set <- proportion[rownames(proportion)%in%gene_list_df$Gene,]
    mean_proportion <- apply(proportion_gene_set,2,mean)
    
    #Permutations
    mean_proportion_bootstrap_df <- matrix(ncol=ncol(proportion_gene_set),nrow=input$perm)
    for (i in 1:input$perm){
      mean_proportion_bootstrap_df[i,] <- apply(proportion[sample(nrow(proportion),nrow(proportion_gene_set),replace=F),],2,mean)
    }

    #Get Pvalue
    
    pvalues <- vector("numeric", ncol(proportion_gene_set))
    for (i in 1:ncol(mean_proportion_bootstrap_df)){
      number_null_more_extreme <- length(which(mean_proportion_bootstrap_df[,i] >= mean_proportion[i]))
      pvalues[i] <- (number_null_more_extreme+1)/(nrow(mean_proportion_bootstrap_df)+1)
    }
    names(pvalues) <- colnames(proportion_gene_set)
    
    #Get Z-score
    
    sd_boot <- apply(mean_proportion_bootstrap_df,2,sd)
    mean_boot <- apply(mean_proportion_bootstrap_df,2,mean)
    z_scores <- (mean_proportion-mean_boot)/sd_boot
    
    results <- as.data.frame(t(rbind(pvalues,z_scores)))
    results <- cbind(cell_type=gsub("\\."," ",rownames(results)),results)
    rownames(results) <- NULL
    results <- arrange(results,pvalues,-z_scores)
    results <- mutate(results,cell_type=factor(cell_type,levels=rev(cell_type)))
    return(results)
  })
  
  #Plot EWCE Cell type enrichment plot
  output$EWCE_plot <- renderPlot({
    req(ewce_results())
    p3 <- ggplot(ewce_results(),aes(cell_type,-log10(pvalues),fill=cell_type)) + geom_bar(stat="identity") + coord_flip() +theme_bw()+ theme(legend.position = "none",text = element_text(size=18)) 
    p3 <- p3 + ylab(expression('-log'[10]*'(pvalue)')) + xlab("") + geom_hline(yintercept = -log10(0.05/nrow(ewce_results())))
    p3
  })
  
  #Save EWCE plot
  output$save_plot_ewce <- downloadHandler(
    filename = function() {
      if(ncol(gene_list_id_matched())>1){
        paste(input$file1,".",input$choose_col,".",min(input$top_n),"-",max(input$top_n),".perm.",input$perm,".EWCE", ".pdf", sep="")
      } else{
        paste(input$file1,"perm",input$perm, "EWCE","pdf", sep=".")
      }
    },
    content = function(file) {
      p3 <- ggplot(ewce_results(),aes(cell_type,-log10(pvalues),fill=cell_type)) + geom_bar(stat="identity") + coord_flip() +theme_bw()+ theme(legend.position = "none",text = element_text(size=18)) 
      p3 <- p3 + ylab(expression('-log'[10]*'(pvalue)')) + xlab("") + geom_hline(yintercept = -log10(0.05/nrow(ewce_results())))
      ggsave(file, plot = p3, device = "pdf")
    }
  )
  
  #Save EWCE table
  output$save_table_ewce <- downloadHandler(
    filename = function() {
      if(ncol(gene_list_id_matched())>1){
        paste(input$file1,".",input$choose_col,".",min(input$top_n),"-",max(input$top_n),".perm.",input$perm,".EWCE", ".txt", sep="")
      } else{
        paste(input$file1,"perm",input$perm,"EWCE", "txt", sep=".")
      }
    },
    content = function(file) {
      write_tsv(ewce_results(), file)
    }
  )
    
})