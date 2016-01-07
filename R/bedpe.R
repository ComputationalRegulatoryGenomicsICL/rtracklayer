# BEDPE Import/Export

# Classes

setClass("BEDPEFile", contains = "RTLFile")

BEDPEFile <- function(resource) {
  new("BEDPEFile", resource = resource)
}

# Export

setGeneric("export.bedpe",
           function(object, con, ...) standardGeneric("export.bedpe"))

setMethod("export.bedpe", "ANY",
          function(object, con, ...) {
            export(object, con, "bedpe", ...)
          })

setMethod("export", c("ANY", "BEDPEFile"),
          function(object, con, format, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format) # needs to be better
            cl <- class(object)
            track <- try(as(object, "GInteractions"), silent = TRUE)
            if (class(track) == "try-error") {
              stop("cannot export object of class '", cl, "': ", track)
          }
            export(track, con, ...)
          })

setMethod("export", c("GInteractions", "BEDPEFile"),
          function(object, con, format, append = FALSE, ignore.strand = FALSE)
          {
            print("inside GI export method")
            if (!missing(format))
              checkArgFormat(con, format)
            object_anchors = anchors(object)
            df <- data.frame(
                             seqnames(object_anchors$first),
                             start(object_anchors$first) - 1L,
                             end(object_anchors$first),
                             seqnames(object_anchors$second),
                             start(object_anchors$second) - 1L,
                             end(object_anchors$second)
                            )
            score <- object$score
            if (!is.null(score)) {
              if (!is.numeric(score) || any(is.na(score)))
                stop("Scores must be non-NA numeric values")
            } else {
                score <- 0
            }
            name <- object$name
            if (is.null(name))
              name <- names(object)
            if (is.null(name))
              name <- rep(NA, length(object))
            df$name <- name
            df$score <- score
            df$strand1 <- strand(object_anchors$first)
            df$strand2 <- strand(object_anchors$second)
            df <- cbind(df, mcols(object))
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            file <- con
            con <- connection(con, if (append) "a" else "w")
            on.exit(release(con))
            write.table(df, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE, na = ".")
            release(con)
          })

# Import

setGeneric("import.bedpe", function(con, ...) standardGeneric("import.bedpe"))

setMethod("import", "BEDPEFile",
          function(con, format, text,
                   genome = NA, colnames = NULL,
                   seqinfo = NULL, extraCols = character())
          {
            if (!missing(format))
              checkArgFormat(con, format)
            file <- con
            con <- queryForConnection(con, which)
            if (attr(con, "usedWhich"))
              which <- NULL
            if (is.null(seqinfo))
              seqinfo <- attr(con, "seqinfo")
            ## check for a track line
            bedpeNames <- c("chrom1", "start1", "end1",
                            "chrom2", "start2", "end2",
                            "name", "score", "strand1", "strand2")
            bedpeClasses <- c("character", "integr", "integer",
                            "character", "integer", "integer",
                            "character", "numeric", "character", "character")
            normArgColnames <- function(validNames) {
              if (is.null(colnames))
                colnames <- validNames
              else {
                colnames <- unique(c(head(bedNames, 3), as.character(colnames)))
                missingCols <- setdiff(colnames, validNames)
                if (length(missingCols))
                  stop("Requested column(s) ",
                       paste("'", missingCols, "'", sep = "", collapse = ", "),
                       " are not valid columns or were not found in the file")
              }
              colnames
            }
            ## read a single line to get ncols up-front,
            ## and thus specify all col classes
            ## FIXME: reading in 'as.is' to save memory,
            line <- readLines(con, 1, warn=FALSE)
            ## UCSC seems to use '#' at beginning to indicate comment.
            while(length(line) &&
                  (!nzchar(line) || substring(line, 1, 1) == "#"))
            {
              line <- readLines(con, 1, warn=FALSE)
            }
            if (length(line)) {
              `tail<-` <- function(x, n, value)
                if (n != 0) c(head(x, -n), value) else x
              pushBack(line, con)
              colsInFile <- seq_len(length(strsplit(line, "[\t ]")[[1]]))
              presentNames <- bedNames[colsInFile]
              tail(presentNames, length(extraCols)) <- names(extraCols)
              bedNames <- presentNames
              presentClasses <- bedClasses[colsInFile]
              tail(presentClasses, length(extraCols)) <- unname(extraCols)
              colnames <- normArgColnames(presentNames)
              bedClasses <- ifelse(presentNames %in% colnames,
                                   presentClasses, "NULL")
              bed <- DataFrame(read.table(con, colClasses = bedClasses,
                                          as.is = TRUE, na.strings = ".",
                                          comment.char = ""))
            } else {
              if (is.null(colnames))
                colnames <- character()
              else colnames <- normArgColnames(bedNames)
              keepCols <- bedNames %in% colnames
              bed <- DataFrame(as.list(sapply(bedClasses[keepCols], vector)))
            }
            colnames(bed) <- bedNames[bedNames %in% colnames]
            bed <- bed[substring(bed$chrom, 1, 1) != "#",]
            # manually handle genome/seqinfo since GenomicData contsructor unavailable
            if (!is.null(seqinfo)) {
                if (is.na(genome))
                    genome <- singleGenome(genome(seqinfo))
                else if (!all(genome == genome(seqinfo), na.rm=TRUE))
                    stop("'genome' ", genome, "' does not match that in 'seqinfo'")
            }
            if (is.null(seqinfo) && !is.na(genome))
                seqinfo <- seqinfoForGenome(genome)
            GInteractions(
              GenomicRanges(bed$chrom1, IRanges(bed$start1, bed$end1), bed$strand1),
              GenomicRanges(bed$chrom2, IRanges(bed$start2, bed$end2), bed$strand2),
              bed[-(c(1:6, 9:10))])
            # TODO Add seqinfo for GInteractions
          })

