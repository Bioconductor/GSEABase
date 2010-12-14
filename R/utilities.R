## Placeholder 'till something appropriate decided
.uniqueIdentifier <- function()
{
    paste(Sys.info()['nodename'], Sys.getpid(), Sys.time(), sep=":")
}

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

.stopf <- function(...)
    stop(simpleError(sprintf(...),
                     call=match.call(call=sys.call(sys.parent(1)))))

.warningf <- function(...)
    warning(simpleWarning(sprintf(...),
                          call=match.call(call=sys.call(sys.parent(1)))))
