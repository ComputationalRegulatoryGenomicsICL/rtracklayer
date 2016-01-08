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
            if (ignore.strand) {
                strand1 <- NA # strand cannot be null in bedpe
                strand2 <- NA
            } else {
                strand1 <- strand(object_anchors$first)
                strand2 <- strand(object_anchors$second)
                strand1[strand1 == "*"] <- NA
                strand2[strand2 == "*"] <- NA
            }
            df$strand1 = strand1
            df$strand2 = strand2
            mcol_names = setdiff(names(mcols(object)), c("name", "score"))
            df <- cbind(df, mcols(object)[,mcol_names])
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
                   genome = NA, seqinfo = NULL, 
                   extraCols = character())
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
            bedpeClasses <- c("character", "integer", "integer",
                            "character", "integer", "integer",
                            "character", "numeric", "character", "character")
            if (length(extraCols) > 0) {
                bedpeNames = c(bedpeNames, names(extraCols))
                bedpeClasses = c(bedpeClasses, unname(extraCols))
            }
            bedpe <- DataFrame(read.table(con, colClasses = bedpeClasses,
                                        as.is = TRUE, na.strings = "."))
            colnames(bedpe)[1:length(bedpeNames)] <- bedpeNames # prevent recycling
            # manually handle genome/seqinfo/strand since GenomicData contsructor unavailable
            bedpe$strand1[is.na(bedpe$strand1)] <- "*"
            bedpe$strand2[is.na(bedpe$strand2)] <- "*"
            if (!is.null(seqinfo)) {
                if (is.na(genome))
                    genome <- singleGenome(genome(seqinfo))
                else if (!all(genome == genome(seqinfo), na.rm=TRUE))
                    stop("'genome' ", genome, "' does not match that in 'seqinfo'")
            }
            if (is.null(seqinfo) && !is.na(genome))
                seqinfo <- seqinfoForGenome(genome)
            out = GInteractions(
              GRanges(bedpe$chrom1, IRanges(bedpe$start1, bedpe$end1), bedpe$strand1, seqinfo=seqinfo),
              GRanges(bedpe$chrom2, IRanges(bedpe$start2, bedpe$end2), bedpe$strand2, seqinfo=seqinfo)
            )
            mcols(out) = bedpe[-(c(1:6, 9:10))]
            out
          })

