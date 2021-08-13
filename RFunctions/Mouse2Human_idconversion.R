#' Mouse2Human
#'
#' @param MouseGenes list of mouse genes
#' @description Credits mainly to Dave Tang's blog (https://www.biostars.org/p/147484/)
#' @importFrom biomaRt useMart getLDS
#' @return table of gene identifiers
#' @export
#'
Mouse2Human <- function(MouseGenes){
	human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	genesMousetoHuman = biomaRt::getLDS(attributes = c("ensembl_gene_id","mgi_symbol", "entrezgene_id"),
														 filters = "mgi_symbol",
														 values = MouseGenes ,
														 mart = mouse,
														 attributesL = c("ensembl_gene_id", "hgnc_symbol"),
														 martL = human,
														 uniqueRows = TRUE)

	colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Entrez.Gene_ID", "Human.Gene_ID", "HGNC")
	return(genesMousetoHuman)
}

