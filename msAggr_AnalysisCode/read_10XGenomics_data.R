read_10XGenomics_data<-function(sample.list){
	data <- sapply(sample.list, FUN = function(x) {x = dirname(list.files(path = x, pattern = 'barcodes.tsv', full.names = TRUE, recursive = TRUE, no.. = TRUE))})
	names(data) <- sample.list
	return(data)
}
