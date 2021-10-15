#' metadata_overlap_stas
#'
#' @description Generates table describing percentage of one metadata col's ids in another metadata col's ids. Note that is currently hardcoded to require numeric data
#'
#' @param seurat.object object of class Seurat
#' @param target.col.pattern metadata column pattern to match against
#' @param xlsx_filename if not NULL, generate an excel file of the enrichment tables with this name
#' @param test.col.pattern metadata column pattern to test for enrichment in coluns matching target.col.pattern
#' @param print_to_screen whether or not to print enrichment dataframes to screen
#' @param xlsx.sheet.count.start if exporting xlsx, sheet number to start counting at. Ignored if xlsx_filename = NULL
#'
#' @return
#' @export
#'
#' @examples
metadata_overlap_stats <- function(seurat.object, target.col.pattern, test.col.pattern, print_to_screen = TRUE, xlsx_filename = NULL, xlsx.sheet.count.start = 1){
	xlsx.sheet.count <- xlsx.sheet.count.start
	for(test.col in colnames(seurat.object@meta.data)){
		if(grepl(pattern = test.col.pattern, x = test.col) == TRUE){
			test.col.id <- levels(seurat.object@meta.data[, test.col])
			for(target.col in colnames(seurat.object@meta.data)){
				if(grepl(pattern = target.col.pattern, x = target.col) == TRUE){
					print(paste("Comparing", test.col, "and", target.col))
					target.col.id <- levels(seurat.object@meta.data[, target.col])
					enrich.df <- data.frame(matrix(nrow = length(target.col.id), ncol = length(test.col.id)))
					rownames(enrich.df) <- target.col.id
					colnames(enrich.df) <- test.col.id
					seurat.meta.df <- seurat.object@meta.data
					for(col.id in colnames(enrich.df)){
						tot.test.cells <- nrow(seurat.meta.df[seurat.meta.df[[test.col]] == col.id,])
						for(row.id in rownames(enrich.df)){
							num.x <- nrow(seurat.meta.df[seurat.meta.df[[target.col]] == row.id & seurat.meta.df[[test.col]] == col.id,])
							pct.x <- as.integer((num.x / tot.test.cells)*100)
							enrich.df[row.id, col.id] <- pct.x
						}
					}
					# if(!sum(rowSums(enrich.df)) == nrow(seurat.object@meta.data)){print("error! Not all cells have been assigned to a category!")}
					colnames(enrich.df) <- sapply(colnames(enrich.df), function(x) paste0(test.col, "_", x))
					rownames(enrich.df) <- sapply(rownames(enrich.df), function(x) paste0(target.col, "_", x))
					if(print_to_screen == TRUE){
						print(enrich.df)
					}					
					if(!is.null(xlsx_filename)){
						sheet.name <- paste0("sheet", xlsx.sheet.count)
						xlsx::write.xlsx(enrich.df, sheetName = sheet.name, file = xlsx_filename, append = TRUE)}
					xlsx.sheet.count <- xlsx.sheet.count + 1
				}
			}
		}
	}
	
}
