\name{RangedData-methods}
\docType{methods}
\alias{RangedData-methods}

% constructor
\alias{GenomicData}

% accessors
\alias{chrom}
\alias{chrom,RangedData-method}
\alias{chrom,GRanges-method}
\alias{chrom<-}
\alias{chrom<-,RangedData-method}
\alias{chrom<-,GRanges-method}
\alias{score,ANY-method}
\alias{score,GenomicRangesORGRangesList-method}

\title{Data on a Genome}

\description{The \code{rtracklayer} package adds
convenience methods on top of \code{RangedData} and \code{GenomicRanges} to
manipulate data on genomic ranges. For \code{RangedData} the spaces are
now called chromosomes (but could still refer to some other type of sequence).
Similarly the universe refers to the genome.}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{RangedData} or \code{GenomicRanges} object.

  \describe{
    \item{}{
      \code{chrom(x), chrom(x) <- value}: Gets or
      sets the chromosome names for \code{x}. The length of \code{value}
      should equal the length of \code{x}. This is an alias for
      \code{\link[IRanges:RangedData-class]{names}(x)}.
    }
    \item{}{
      \code{score(x)}: Gets the \dQuote{score} column from the element
      metadata of a \code{GenomicRanges} or \code{GRangesList}. Many
      track formats have a score column, so this is often used during
      export. The IRanges package defines a method for
      \code{RangedData}. The \code{ANY} fallback for this method simply
      returns \code{NULL}.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicData(ranges, ..., strand = NULL, chrom = NULL,
        genome = NULL)}: Constructs a \code{GRanges} instance with
      the given \code{ranges} and variables in \code{...} (see the
      \code{\link[GenomicRanges]{GRanges}} constructor).

      If non-\code{NULL}, the \code{strand} argument specifies the strand
      of each range. It should be a character vector or factor of length
      equal to that of \code{ranges}. All values should be either \code{-},
      \code{+}, or \code{*}. To get the levels for \code{strand}, call
      \code{levels(strand())}.

      \code{chrom} argument is analogous to \code{seqnames} in the
      \code{GRanges} and \code{space} in \code{RangedData} constructors.

      The \code{genome} argument should be a scalar string and is treated
      as the \code{RangedData} universe. See the examples.
      
      If \code{ranges} is not a \code{Ranges} object, this function
      calls \code{as(ranges, "RangedData")} and returns the result if
      successful. As a special case, the \dQuote{chrom} column in a
      \code{data.frame}-like object is renamed to \dQuote{space}, for
      convenience. Thus, one could pass a \code{data.frame} with columns
      \dQuote{start}, \dQuote{end} and, optionally, \dQuote{chrom}.
    }
  }
}

\author{ Michael Lawrence and Patrick Aboyoun }

\examples{
  range1 <- IRanges::RangesList(chr1 = IRanges::IRanges(c(1,2,3), c(5,2,8)))

  ## just ranges ##
  ## GRanges instance
  gr <- GenomicData(range1)

  ## with a genome (universe) ##
  ## GRanges instance
  gr <- GenomicData(range1, genome = "hg18")
  genome(gr) ## "hg18"

  ## with some data ##
  filter <- c(1L, 0L, 1L)
  score <- c(10L, 2L, NA)
  strand <- factor(c("+", NA, "-"), levels = levels(strand()))
  ## GRanges instance
  gr <- GenomicData(range1[[1]], score, chrom = "chr1", genome = "hg18")
  mcols(gr)[["score"]]
  strand(gr) ## all '*'
  gr <- GenomicData(range1[[1]], score, filt = filter, strand = strand,
                    chrom = "chr1")
  mcols(gr)[["filt"]]
  strand(gr) ## equal to 'strand'
  ## coercion from data.frame ##
  df <- as.data.frame(gr)
}
\keyword{classes}
\keyword{methods}
