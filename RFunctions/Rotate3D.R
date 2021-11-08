#' Rotate3D
#'
#' @param embedded.matrix 
#' @param output.dir 
#' @param color.palette 
#' @param window.rect 
#' @param image.zoom 
#'
#' @return
#' @export
#'
#' @examples
rotate_3d <- function(embedded.matrix, output.dir = "Rotated_3D_Animation", color.palette = NULL, window.rect = c(122, 44, 950, 929), image.zoom = 0.6){	
	my_object3d<-data.matrix(embedded.matrix)
	scatterplot3d(x=my_object3d[,1], y=my_object3d[,2], z=my_object3d[,3], color=c("red"))
	if(!is.null(color.palette)){
	plot3d(x=my_object3d[,1], y=my_object3d[,2], z=my_object3d[,3], 
				 col = color.palette, 
				 box=FALSE, axes=FALSE, xlab = "", ylab = "", zlab = "")
	} else {
		plot3d(x=my_object3d[,1], y=my_object3d[,2], z=my_object3d[,3], 
					 box=FALSE, axes=FALSE, xlab = "", ylab = "", zlab = "")
	}
	par3d("windowRect" = window.rect)
	rgl.viewpoint(zoom=image.zoom)
	animation_dir<-paste0(output.dir)
	try(dir.create(animation_dir))
	for (i in 1:90) {
		view3d(userMatrix=rotationMatrix(2*pi * i/90, 1, -1, -1), zoom=image.zoom)
		rgl.snapshot(filename=paste0(animation_dir,"/frame-",
																sprintf("%03d", i), ".png"))
	}
	
}
