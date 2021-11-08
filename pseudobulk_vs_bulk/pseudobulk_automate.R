## Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(xlsx)

## Set global variables
projectName <- "pseudobulk_allPops"
cell.list <- c("LSKm2", "CMPm2", "CMPm", "MEPm", "GMPm")
bulk.list <- c("LSK", "CMP", "CFUE", "CFUMK", "ERY", "GMP", "MK", "MEP")
bulk.folders <- c("TotalSeq", "ScriptSeq")
count.progress <- 0

## Set function options
plot.graphs <- FALSE
export.graphs <- FALSE
all.comparisons <- TRUE

## Store session info
sessionInfo.filename <- paste0(projectName, "_sessionInfo.txt")
sink(sessionInfo.filename)
sessionInfo()
sink()

## Create empty spearman coefficient dataframes
spearman.scriptseq <- data.frame(matrix(NA, nrow = length(bulk.list), ncol = length(cell.list)))
rownames(spearman.scriptseq) <- bulk.list
colnames(spearman.scriptseq) <- cell.list
spearman.totalseq <- spearman.scriptseq

if(plot.graphs == FALSE){
	print("Export graphs set to FALSE")
}
if(export.graphs == FALSE){
	print("Export graphs set to FALSE")
}

for (cell.pop in cell.list){
	# Load scRNASeq object
	seurat.object <- readRDS(paste0("~/Desktop/10XGenomicsData/cellRanger/", cell.pop, "_raw.RDS"))
	# Calculate cpm 
	sc.mtx <- data.frame(rowSums(seurat.object@assays$RNA@counts))
	sc.mtx <- sc.mtx/colSums(sc.mtx)*1000000
	
	if(all.comparisons == TRUE){
		bulk.pops <- bulk.list
	}else{
		bulk.pops <- substr(cell.pop, 1, 3)
	}
	for(comparison.pop in bulk.pops){
		# Automate through bRNASeq data
		
		for(sequencing.folder in bulk.folders){
			bulk.subfolders <- c(paste0("~/Desktop/10XGenomicsData/bRNASeq/", sequencing.folder, "/"))
			
			for(bulk.path in bulk.subfolders){
				bulk.path <- bulk.subfolders[1]
				file_list <- list.files(path = bulk.path, pattern = comparison.pop)
				bulk.pop <- read.table(file = paste0(bulk.path, file_list[1]), header = TRUE)[1]
				
				for(file in file_list){
					bulk.pop <- merge(bulk.pop, read.table(file = paste0(bulk.path, file), header = TRUE, sep = "\t")[, c(1,5)], by = "gene_id")
				}
				
				#calculate average tpm
				bulk.pop[comparison.pop] <- rowMeans(x = bulk.pop[, c(2, 3)])
				
				# clean and check df
				bulk.pop <- bulk.pop[, -c(2, 3)]
				
				
				#calc log_rho's
				bulk.rho <- log((bulk.pop[comparison.pop] + 1), 2) - log(colSums((bulk.pop[comparison.pop])+ 1), 2)
				rownames(bulk.rho) <- bulk.pop$gene_id
				colnames(bulk.rho) <- "bulk"
				sc.rho <- log((sc.mtx + 1), 2) - log(colSums((sc.mtx+1)), 2)
				colnames(sc.rho) <- "sc"
				
				# check all entries are numeric
				if(is.numeric(bulk.rho[,1]) == FALSE){
					print("WARNING: bulk.rho datafame contains non-numeric values")
				}
				if(is.numeric(sc.rho[,1]) == FALSE){
					print("WARNING: sc.rho datafame contains non-numeric values")
				}
				
				
				# drop NAs
				concordance.df <- merge(bulk.rho, sc.rho, by = "row.names", all = TRUE)
				concordance.df <- na.omit(concordance.df)
				rownames(concordance.df) <- concordance.df$Row.names
				concordance.df <- concordance.df[, c(2, 3)]
				
				# spearman plot
				spearman.coeff <- cor(concordance.df, method = "spearman")[1, 2]
				# spearman.coeff <- round(spearman.coeff, 3)
				if(sequencing.folder == "TotalSeq"){
					spearman.totalseq[comparison.pop, cell.pop] <- spearman.coeff
				}
				if(sequencing.folder == "ScriptSeq"){
					spearman.scriptseq[comparison.pop, cell.pop] <- spearman.coeff
				}
				
				
				if(plot.graphs == TRUE){
					cor.plot <- ggplot(data = concordance.df, aes(x = bulk, y = sc)) +
						geom_point()  +
						ggtitle(paste0(cell.pop, "_vs_", comparison.pop, "-", sequencing.folder, ": spearman coeff = ", round(spearman.coeff, digits = 3)))
					plot(cor.plot)
				}
				
				if(export.graphs == TRUE){
					png(filename = paste0(cell.pop, "_vs_", comparison.pop, "-", sequencing.folder,"-spearmanPlot.png"), height = 480, width = 800)
					plot(cor.plot)
					dev.off()
				}
			}
			count.progress <- count.progress + 1
			print(count.progress)
			
		}
	}
}

spearman.scriptseq
spearman.totalseq
spearman.df

write.xlsx(round(spearman.scriptseq, 3), file = "pseudobulk_allpops_spearmanCoeff.xlsx", sheetName = "ScriptSeq", col.names = TRUE, row.names = TRUE)
write.xlsx(round(spearman.totalseq, 3), file = "pseudobulk_allpops_spearmanCoeff.xlsx", sheetName = "TotalSeq", col.names = TRUE, row.names = TRUE, append = TRUE)





