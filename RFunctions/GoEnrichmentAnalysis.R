# #####NOTE: This program requires internet access

#' GO_enrichment_function
#'
#' @param my_folder folder containing files with gene lists in ensembl or symbol format
#' @param number_of_terms how many terms to return (use 'Inf' to return everything)
#' @param list.species currently only mm10 is supported
#' @param ontology.list ontology types to return (bioprocess, molecular function, or cell compartment)
#' @param file.pattern grep suffix pattern describing files in folder that contain gene lists
#' @param outfile.suffix suffix to add to base file name
#' @param is.header Logical; does gene file have a header
#' @param id.type gene ID type (either ensmembl or symbol)
#'
#' @importFrom limma goana topGO
#' @importFrom org.Mm.eg.db org.Mm.egENSEMBL org.Mm.egALIAS2EG
#' @return
#' @export
#'
GO_enrichment_function <- function(my_folder, number_of_terms = Inf, list.species = 'Mm', ontology.list = c('BP', 'MF', 'CC'), file.pattern = '.txt$', outfile.suffix = "-GOAll.txt", id.type = 'ensembl', is.header = FALSE){

	# setwd(my_folder)
	file.list <- list.files(path = my_folder, pattern = file.pattern, full.names = T)

	for(i in file.list){
		outfile = gsub(file.pattern, "", basename(i))
		print(paste("creating", outfile))
		if(is.header == TRUE){
			genelist<-read.delim(file = i, header = TRUE)
		}
		if(is.header == FALSE){
			genelist<-read.delim(file = i, header = FALSE)
		}
		
		if(tolower(id.type) == 'ensembl'){
			# No geneID conversion required -------------------------------------------
			table_list_entrezIDs <- unlist(mget(x=as.character(genelist[,1]), envir=org.Mm.eg.db::org.Mm.egENSEMBL, ifnotfound=NA))
		}	else if(tolower(id.type) == 'symbol'){
			# Symbol to Entrez conversion --------------------------------------------
			table_list_entrezIDs <- unlist(mget(x=as.character(genelist[,1]), envir=org.Mm.eg.db::org.Mm.egALIAS2EG, ifnotfound=NA))
		} else {
			stop ("ID not recognized")
		}

		table_list_goana <- limma::goana(table_list_entrezIDs, species=list.species )
		table_list_GO <- limma::topGO(table_list_goana, ontology=ontology.list, number = number_of_terms)
		try(write.table(table_list_GO, file=paste(my_folder, outfile, outfile.suffix, sep = ""), quote = F, row.names = T, sep = "\t"))
		try(remove(genelist, table_list_entrezIDs, table_list_goana, table_list_GO))
	}
}
