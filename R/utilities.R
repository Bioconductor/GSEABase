## simplified unique for vectors, preserving attributes
.unique <- function(x, y) c(x, y[!(y %in% x)])

.glue <- function(x, y, op)
    paste("(", x, " ", op, " ", y, ")", sep="")

