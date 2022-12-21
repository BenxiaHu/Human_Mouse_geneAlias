rm(list=ls())
library(DT)

createLink <- function(base,val) {
  sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>',base,val) ##target="_blank" 
}

getGeneName <- function(ref1,ref2,ref3,flag){
  eg2symbol <- toTable(ref1)
  eg2name <- toTable(ref2)
  eg2alias <- toTable(ref3)
  
  eg2alis_list <- lapply(split(eg2alias,eg2alias$gene_id),function(x){paste0(x[,2],collapse = ";")})
  GeneList <- mappedLkeys(ref1)
  
  if( GeneList[1] %in% eg2symbol$symbol ){
    symbols <- GeneList
    geneIds <- eg2symbol[match(symbols,eg2symbol$symbol),'gene_id']
  }else{
    geneIds <- GeneList
    symbols <- eg2symbol[match(geneIds,eg2symbol$gene_id),'symbol']
  }
  geneNames <- eg2name[match(geneIds,eg2name$gene_id),'gene_name']
  geneAlias <- sapply(geneIds,function(x){ifelse(is.null(eg2alis_list[[x]]),"no_alias",eg2alis_list[[x]])})
  
  gene_info <- data.frame(symbols=symbols,
                          geneIds=createLink(paste0("http://www.ncbi.nlm.nih.gov/gene/",geneIds),geneIds),
                          geneNames=geneNames,
                          geneAlias=geneAlias,
                          stringsAsFactors = F
  )
  colnames(gene_info) <- paste0(flag,"_",c("symbols","geneIds","geneNames","geneAlias"))
  return(gene_info)
}

library(org.Mm.eg.db)
Mousegene <- getGeneName(org.Mm.egSYMBOL,org.Mm.egGENENAME,org.Mm.egALIAS2EG,"Mouse")
detach(package:org.Mm.eg.db,unload=TRUE)
Mousegene$symbol <- toupper(Mousegene$Mouse_symbols)

library(org.Hs.eg.db)
Humangene <- getGeneName(org.Hs.egSYMBOL,org.Hs.egGENENAME,org.Hs.egALIAS2EG,"Human")
detach(package:org.Hs.eg.db,unload=TRUE)

######
library(tidyverse)

result <- full_join(Humangene,Mousegene, by = c("Human_symbols" = "symbol"), keep = F)

#library("xtable") 
#print(xtable(gene_info), type="html",include.rownames=F, file='all_gene.anno',sanitize.text.function = force)

y <- DT::datatable(result,escape = F,rownames=F,
                   options = list(pageLength = 20 ))

file <- 'human_mouse_all_gene_bioconductor.html'
DT::saveWidget(y,file) 


