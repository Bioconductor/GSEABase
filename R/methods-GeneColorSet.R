.constructors("GeneColorSet",
              required=c("setName", "setIdentifier", "phenotype"))

setMethod("initialize",
          signature=signature(.Object="GeneColorSet"),
          function(.Object, .Template=.Object, ...,
                   ## additional args
                   genes=.Template@genes,
                   phenotype=.Template@phenotype,
                   geneColor=factor(character(length(genes))),
                   phenotypeColor=factor(character(length(genes)))) {
              callNextMethod(.Object, .Template, ...,
                             genes=genes,
                             phenotype=mkScalar(phenotype),
                             geneColor=geneColor,
                             phenotypeColor=phenotypeColor)
          })

setMethod("GeneColorSet",
          signature=signature(type="GeneSet"),
          function(type, ...) {
              new("GeneColorSet", type, ...)
          })

setAs("GeneSet", "GeneColorSet",
      function(from)
      stop("use 'GeneColorSet(<GeneSet>, phenotype=\"phenotype\", ...)'"))

.GETTERS_GeneColorSet <- c("phenotype", "geneColor", "phenotypeColor")

.getters("GeneColorSet", .GETTERS_GeneColorSet)

.SETTERS_GeneColorSet <- .GETTERS_GeneColorSet

.setters("GeneColorSet", .SETTERS_GeneColorSet)

setReplaceMethod("phenotype",
                 signature(object="GeneColorSet",
                           value="character"),
                 function(object, value) {
                     slot(object, "phenotype") <- mkScalar(value)
                     object
                 })

setMethod("coloring",
          signature=signature(object="GeneColorSet"),
          function(object, ...) {
              data.frame(geneColor=geneColor(object),
                         phenotypeColor=phenotypeColor(object),
                         row.names=genes(object))
          })

setReplaceMethod("coloring",
                 signature=signature(
                   object="GeneColorSet",
                   value="data.frame"),
                 function(object, ..., value) {
                     ogenes <- genes(object)
                     if (!all(row.names(value) %in% ogenes))
                         stop("'data.frame' row.names must all be gene symbols")
                     if (nrow(value) != length(ogenes))
                         stop("'data.frame' must define colors for all genes")
                     if (length(colnames(value)) !=2 ||
                                !all(c("geneColor", "phenotypeColor") %in%
                                     colnames(value)))
                         stop("'data.frame' must only 'geneColor' and 'phenotypeColor' columns")
                     if (!all("factor" %in% sapply(value, class)))
                         stop("'data.frame columns must be of class 'factor'")
                     new(class(object), object,
                         genes=genes(object),
                         geneColor=value[ogenes, "geneColor"],
                         phenotypeColor=value[ogenes, "phenotypeColor"],
                         setName=setName(object),
                         setIdentifier=setIdentifier(object),
                         phenotype=phenotype(object))
                 })

## Logic operations

.checkGeneColorSetLogicTypes <- function(x, y, functionName) {
    .checkGeneSetLogicTypes(x, y, functionName)
    if (phenotype(x) != phenotype(y))
        warning(functionName, ": ",
                "'phenotype' differs; creating synthetic phenotype")
    else {
        if (any(levels(geneColor(x)) != levels(geneColor(y))) ||
            any(levels(phenotypeColor(x)) !=
                levels(phenotypeColor(y))))
            warning(functionName, ": ",
                    "'levels' of gene- or phenotypeColor differ between identical phenotpyes")
    }
}

.geneColorSetIntersect <- function(x, y) {
    color <- function(x, y, lbl) {
        if (!phenotypesIdentical ||
            any(as.character(x) != as.character(y)))
            factor(.glue(as.character(x), as.character(y), ", "))
        else x
    }
    .checkGeneColorSetLogicTypes(x, y, "'&' or 'intersect'")
    vx <- as.vector(genes(x))
    vy <- as.vector(genes(y))
    idx <- match(vy, vx, 0)             # x index
    idy <- match(vx[idx], vy, 0)        # y index
    genes <- vy[idy]
    phenotype <- phenotype(x)
    phenotypesIdentical <- phenotype == phenotype(y)
    if (!phenotypesIdentical)
        phenotype <- .glue(phenotype, phenotype(y), ", ")
    gc <- color(geneColor(x)[idx], geneColor(y)[idy], "geneColor")
    pc <- color(phenotypeColor(x)[idx], phenotypeColor(y)[idy],
                "phenotypeColor")
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName = .glue(setName(x), setName(y), " & "),
        urls = .unique(urls(x), urls(y)),
        phenotype = phenotype,
        genes = genes, geneColor = gc, phenotypeColor = pc)      
}

.geneColorSetUnion <- function(x, y) {
    .checkGeneColorSetLogicTypes(x, y, "'|' or 'union'")
    idy <- which(!(genes(y) %in% genes(x)))
    genes <- c(genes(x), genes(y)[idy])
    phenotype <- phenotype(x)
    phenotypesIdentical <- phenotype == phenotype(y)
    if (!phenotypesIdentical)
        phenotype <- .glue(phenotype, phenotype(y), ", ")
    gc <- factor(c(as.character(geneColor(x)),
                   as.character(geneColor(y)[idy])),
                 levels = unique(c(
                   levels(geneColor(x)),
                   levels(geneColor(y)))))
    pc <- factor(c(as.character(phenotypeColor(x)),
                   as.character(phenotypeColor(y)[idy])),
                 levels = unique(c(
                   levels(phenotypeColor(x)),
                   levels(phenotypeColor(y)))))
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName = .glue(setName(x), setName(y), " | "),
        urls = .unique(urls(x), urls(y)),
        genes = genes, geneColor = gc, phenotypeColor = pc)      
}

.geneColorSetSetdiff <- function(x, y) {
    .checkGeneColorSetLogicTypes(x, y, "'setdiff'")
    gx <- genes(x)
    gy <- genes(y)
    idx <- 
        if (length(gx) || length(gy))
            match(gx, gy, 0)==0
        else
            TRUE
    gc <- factor(geneColor(x)[idx],
                 levels=levels(geneColor(x)))
    pc <- factor(phenotypeColor(x)[idx],
                 levels=levels(phenotypeColor(x)))
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName = .glue(setName(x), setName(y), " - "),
        urls = .unique(urls(x), urls(y)),
        genes=gx[idx], geneColor = gc, phenotypeColor = pc)
}

setMethod("intersect",
          signature=signature(
            x="GeneColorSet", y="GeneColorSet"),
          .geneColorSetIntersect)

setMethod("union",
          signature=signature(
            x="GeneColorSet", y="GeneColorSet"),
          .geneColorSetUnion)

setMethod("&",
          signature=signature(e1="GeneColorSet", e2="GeneColorSet"),
          function(e1, e2) .geneColorSetIntersect(e1, e2))

setMethod("&",
          signature=signature(e1="GeneColorSet", e2="character"),
          function(e1, e2) {
              idx <- which(genes(e1)==e2)
              new(class(e1), e1,
                  setIdentifier=setIdentifier(e1),
                  setName=.glue(setName(e1), "<character>", " & "),
                  genes=genes(e1)[idx],
                  geneColor=geneColor(e1)[idx],
                  phenotypeColor=phenotypeColor(e1)[idx])
          })

setMethod("|",
          signature=signature(e1="GeneColorSet", e2="GeneColorSet"),
          function(e1, e2) .geneColorSetUnion(e1, e2))

setMethod("|",
          signature=signature(e1="GeneColorSet", e2="character"),
          function(e1, e2) {
              if (!all(e2 %in% e1))
                  stop("named genes not present in ", class(e1))
              e1
          })

setMethod("setdiff",
          signature=signature(
            x="GeneColorSet", y="GeneColorSet"),
          .geneColorSetSetdiff)



        
## other methods

setMethod("show",
          signature=signature(object="GeneColorSet"),
          function(object) {
              callNextMethod()
              cat("phenotype:", phenotype(object), "\n")
              cat("geneColor: ",
                  paste(selectSome(as.character(geneColor(object)),
                                   maxToShow=4),
                        collapse=", "),
                  "\n  levels: ", paste(levels(geneColor(object)),
                                        collapse=", "), "\n",
                  "phenotypeColor: ",
                  paste(selectSome(as.character(phenotypeColor(object)),
                                   maxToShow=4),
                        collapse=", "),
                  "\n  levels: ", paste(levels(phenotypeColor(object)),
                                        collapse=", "), "\n",
                  sep="")
          })
