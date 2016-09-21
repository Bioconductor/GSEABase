## Placeholder 'till something appropriate decided
.uniqueIdentifier <- local({
    node <- NULL
    pid <- NULL
    uid <- 0L
    function() {
        if (is.null(node)) {
            node <<- Sys.info()['nodename']
            pid <<- Sys.getpid()
        }
        uid <<- uid + 1L
        base::paste(node, pid, date(), uid, sep=":")
    }
})

## simplified unique for vectors, preserving attributes
.unique <- function(x, y) c(x, y[!(y %in% x)])

.glue <- function(x, y, op)
    paste("(", x, op, y, ")", sep="")

.requireQ <- function(pkg)
  suppressWarnings(require(pkg, quietly=TRUE, character.only=TRUE))

.nameAll <- function(x) {
    ## Add names to character vector x.  Elements of x without names get
    ## a name matching the element.
    ##
    if (!is.character(x))
      stop("argument 'x' must be a character vector")
    if (length(names(x)))
      names(x) <- ifelse(nchar(names(x)) == 0, x, names(x))
    else
      names(x) <- x
    x
}

.stopf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    stop(simpleError(msg, call=call))
}

.warningf <- function(...) {
    call <- match.call(call=sys.call(sys.parent(1)))
    msg <- paste(sprintf(...), collapse="\n")
    warning(simpleWarning(msg, call=call))
}
