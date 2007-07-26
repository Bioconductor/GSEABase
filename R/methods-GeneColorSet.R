.constructors_GeneSet("GeneColorSet", required=c("phenotype"))

## See methods-GeneSet for rationale about 'initialize'
.colorConstructor <- function(color, len) {
    if (length(color)==0 && len !=0)
        factor(character(len))
    else
        color
}

setMethod("initialize",
          signature=signature(.Object="GeneColorSet"),
          function(.Object, .Template=.Object, ...,
                   ## additional args
                   geneIds=.Template@geneIds,
                   phenotype=.Template@phenotype,
                   geneColor=.Template@geneColor,
                   phenotypeColor=.Template@phenotypeColor) {
              callNextMethod(.Object, .Template, ...,
                             geneIds=geneIds,
                             phenotype=mkScalar(phenotype),
                             geneColor=.colorConstructor(
                               geneColor,
                               length(geneIds)),
                             phenotypeColor=.colorConstructor(
                               phenotypeColor,
                               length(geneIds)))
          })

setMethod("GeneColorSet",
          signature=signature(type="GeneSet"),
          function(type,
                   phenotype,
                   geneColor=factor(character(length(geneIds(type)))),
                   phenotypeColor=factor(character(length(geneIds(type)))),
                   ..., setIdentifier=.uniqueIdentifier()) {
              new("GeneColorSet", as(type, "GeneColorSet"),
                  setIdentifier=setIdentifier,
                  phenotype=phenotype,
                  geneColor=geneColor, phenotypeColor=phenotypeColor,
                  ...)
          })

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
                         row.names=geneIds(object))
          })

setReplaceMethod("coloring",
                 signature=signature(
                   object="GeneColorSet",
                   value="data.frame"),
                 function(object, ..., value) {
                     ogenes <- geneIds(object)
                     if (!all(row.names(value) %in% ogenes))
                         stop("'data.frame' row.names must all be gene symbols")
                     if (nrow(value) != length(ogenes))
                         stop("'data.frame' must define colors for all geneIds")
                     if (length(colnames(value)) !=2 ||
                                !all(c("geneColor", "phenotypeColor") %in%
                                     colnames(value)))
                         stop("'data.frame' must only 'geneColor' and 'phenotypeColor' columns")
                     if (!all("factor" %in% sapply(value, class)))
                         stop("'data.frame columns must be of class 'factor'")
                     new(class(object), object,
                         geneIds=geneIds(object),
                         geneColor=value[ogenes, "geneColor"],
                         phenotypeColor=value[ogenes, "phenotypeColor"],
                         setName=setName(object),
                         setIdentifier=setIdentifier(object),
                         phenotype=phenotype(object))
                 })

## subset

setMethod("[",
          signature=signature(
            x="GeneColorSet", i="numeric"),
          function(x, i, j, ..., drop=TRUE) {
              if (any(duplicated(i)))
                  stop("duplicate index: ",
                       paste(i[duplicated(i)], collapse=" "))
              geneIds <- geneIds(x)[i]
              if (any(is.na(geneIds)))
                  stop("unmatched index: ",
                       paste(i[is.na(geneIds)], collapse=" "))
              new(class(x), x,
                  geneIds=geneIds,
                  geneColor=factor(
                    as.character(geneColor(x)[i])),
                  phenotypeColor=factor(
                    as.character(phenotypeColor(x)[i])))
          })

setMethod("[",
          signature=signature(
            x="GeneColorSet", i="character"),
          function(x, i, j, ..., drop=TRUE) {
              idx <- pmatch(i, geneIds(x))
              if (any(is.na(idx)))
                  stop(sprintf("unmatched / duplicate geneIds: '%s'",
                               paste(i[is.na(idx)], collapse="', '")))
              new(class(x), x,
                  geneIds=geneIds(x)[idx],
                  geneColor=factor(
                    as.character(geneColor(x))[idx]),
                  phenotypeColor=factor(
                    as.character(phenotypeColor(x)[idx])))
          })

setMethod("[[",
          signature=signature(
            x="GeneColorSet", i="numeric"),
          function(x, i, j, ...) {
              c(geneId=geneIds(x)[[i]],
                geneColor= as.character(geneColor(x)[[i]]),
                phenotypeColor= as.character(phenotypeColor(x)[[i]]))
          })

setMethod("[[",
          signature=signature(
            x="GeneColorSet", i="character"),
          function(x, i, j, ...) {
              idx <- match(i, geneIds(x))
              if (is.na(idx))
                  stop("unmatched gene: ", i)
              ## 'next' method is GeneSet, so want to re-start...
              callGeneric(x, idx, ...)
          })

setMethod("$",
          signature=signature(x="GeneColorSet"),
          function(x, name) {
              i <- pmatch(name, geneIds(x), duplicates.ok=FALSE)
              if (is.na(i))
                  stop("unmatched gene: ", i)
              c(geneId=geneIds(x)[i],
                geneColor=as.character(geneColor(x)[[i]]),
                phenotypeColor=as.character(phenotypeColor(x)[[i]]))
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
            any(levels(x) != levels(y)) ||
            any(as.character(x) != as.character(y)))
            factor(.glue(as.character(x), as.character(y), ", "))
        else x
    }
    .checkGeneColorSetLogicTypes(x, y, "'&' or 'intersect'")
    vx <- as.vector(geneIds(x))
    vy <- as.vector(geneIds(y))
    idx <- match(vy, vx, 0)             # x index
    idy <- match(vx[idx], vy, 0)        # y index
    geneIds <- vy[idy]
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
        geneIds = geneIds, geneColor = gc, phenotypeColor = pc)      
}

.geneColorSetUnion <- function(x, y) {
    .checkGeneColorSetLogicTypes(x, y, "'|' or 'union'")
    idy <- which(!(geneIds(y) %in% geneIds(x)))
    geneIds <- c(geneIds(x), geneIds(y)[idy])
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
        geneIds = geneIds, geneColor = gc, phenotypeColor = pc)      
}

.geneColorSetSetdiff <- function(x, y) {
    .checkGeneColorSetLogicTypes(x, y, "'setdiff'")
    gx <- geneIds(x)
    gy <- geneIds(y)
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
        geneIds=gx[idx], geneColor = gc, phenotypeColor = pc)
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
              idx <- which(geneIds(e1)==e2)
              new(class(e1), e1,
                  setIdentifier=setIdentifier(e1),
                  setName=.glue(setName(e1), "<character>", " & "),
                  geneIds=geneIds(e1)[idx],
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
                  stop("named geneIds not present in ", class(e1))
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
