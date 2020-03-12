flip_strand <- function(x) {
	x <- gsub("A", "V", x)
	x <- gsub("T", "X", x)
	x <- gsub("C", "Y", x)
	x <- gsub("G", "Z", x)
	x <- gsub("V", "T", x)
	x <- gsub("X", "A", x)
	x <- gsub("Y", "G", x)
	x <- gsub("Z", "C", x)
	return(x)   
}
