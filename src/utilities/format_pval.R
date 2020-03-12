# Custom p-value formatting
my_format_pval <- function(p) {
  sapply(p, function(x) {
		if (is.na(x)) {
			return("")
		} else if (x >= 0.06) {
			round(x, digits=2)
		} else if (x >= 0.001) {
			round(x, digits=3)
		} else {
			format(x, digits=1, scientific=TRUE)
		}
  })
}
