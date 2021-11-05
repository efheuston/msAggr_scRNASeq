#' Metadata.Stats
#' @description for each metadata column matching regex X, what's that number and percentage of cells in each?
#'
#' @param seurat.object object of class Seurat
#' @param metadata.pattern grepl pattern to match for metadata columns
#'
#' @return
#' @export
#'
#' @examples
metadata_stats <- function(seurat.object, metadata.pattern){
	tot.cells <- nrow(seurat.object@meta.data)
	for (meta.col in colnames(cmp.object@meta.data)){
		if(grepl(pattern = metadata.pattern, x = meta.col)==TRUE){
			new.clusters <- sort(as.numeric(levels(seurat.object@meta.data[[meta.col]])))
			stats.df <- data.frame(matrix(ncol = 2, nrow = length(new.clusters)))
			colnames(stats.df) <- c("num_cells", "pct_pop")
			rownames(stats.df) <- new.clusters
			meta.df <- seurat.object@meta.data
			for(row.id in rownames(stats.df)){
				num.x <- nrow(meta.df[meta.df[meta.col] == row.id,])
				pct.x <- as.integer(num.x / tot.cells *100)
				# print(pct.x)
				stats.df[row.id, "num_cells"] <- num.x
				stats.df[row.id, "pct_pop"] <- pct.x
			}
			print(stats.df)
		}
	}
}
