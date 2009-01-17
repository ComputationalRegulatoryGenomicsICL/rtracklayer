\name{RangedData-methods}
\docType{methods}
\alias{RangedData-methods}

% constructor
\alias{GenomicData}

% accessors
\alias{score<-}
\alias{score<-,RangedData-method}
\alias{score,RangedData-method}
\alias{chrom}
\alias{chrom,RangedData-method}
\alias{genome}
\alias{genome,RangedData-method}
\alias{genome<-}
\alias{genome<-,RangedData-method}
\alias{strand,RangedData-method}

\title{Data on a Genome} \description{The \code{rtracklayer} package adds
convenience methods on top of \code{RangedData} manipulate data on
genomic ranges. The spaces are now called chromosomes (but could still
refer to some other type of sequence). The universe refers to the genome
and there is formal treatment of an optional \code{strand} variable.}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{RangedData} object.

  \describe{
    \item{}{
      \code{chrom(x)}: Gets the chromosome names (a factor) over the
      ranges in \code{x}.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}: Gets or sets the
      genome (a single string or \code{NULL}) for
      the ranges in \code{x}; simple wrappers around
      \code{\link[IRanges]{universe}} and
      \code{\link[IRanges]{universe<-}}, respectively.
    }
    \item{}{
      \code{\link[BSgenome]{strand}(x)}: Gets the strand factor in
      \code{x}, with the standard 
      levels: \code{-}, \code{+} and \code{*}, referring to the negative,
      positive and either/both strands, respectively. Any strand factor
      stored in the \code{RangedData} should have those levels.
      \code{NA}'s are allowed; the value is all \code{NA} if no strand
      has been specified.
    }
    \item{}{\code{score(x)}: gets the column representing a "score" in
      \code{x}, as a vector. This is the column named
      \code{score}, or, if this does
      not exist, the first column, if it is numeric. Otherwise,
      \code{NULL} is returned.
    }
    \item{}{\code{score(x) <- value}: sets the column named \code{score}
      to \code{value}, which should be a numeric vector of length equal
      to the number of rows.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicData(ranges, ..., strand = NULL, chrom = NULL,
        genome = NULL)}: Constructs a \code{RangedData} instance with
      the given \code{ranges} and variables in \code{...} (see the
      \code{\link[IRanges]{RangedData}} constructor). If
      non-\code{NULL}, \code{strand} specifies the strand of each
      range. It should be a character vector or factor of length equal to
      that of \code{ranges}. All values should be either \code{-}, \code{+},
      \code{*} or \code{NA}. To get these levels, call
      \code{levels(strand())}. \code{chrom} is analogous to
      \code{splitter} in \code{RangedData}; if non-\code{NULL} it
      should be coercible to a factor indicating how the ranges,
      variables and strand should be split up 
      across the chromosomes. The \code{genome} argument should be a
      scalar string and is treated as the \code{RangedData}
      universe. See the examples.
    }
  }
}

\author{ Michael Lawrence }

\examples{
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))

  ## just ranges
  gr <- GenomicData(range1) 

  ## with a genome (universe)
  gr <- GenomicData(range1, genome = "hg18")
  genome(gr) ## "hg18"

  ## with some data
  filter <- c(1L, 0L, 1L)
  score <- c(10L, 2L, NA)
  strand <- factor(c("+", NA, "-"), levels = levels(strand()))
  gr <- GenomicData(range1, score, genome = "hg18")
  gr[["score"]]
  strand(gr) ## all NA
  gr <- GenomicData(range1, score, filt = filter, strand = strand)
  gr[["filt"]]
  strand(gr) ## equal to 'strand'

  range2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  ranges <- c(range1, range2)
  score <- c(score, c(0L, 3L, NA, 22L)) 
  chrom <- paste("chr", rep(c(1,2), c(length(range1), length(range2))), sep="")
  
  gr <- GenomicData(ranges, score, chrom = chrom, genome = "hg18")
  chrom(gr) # equal to 'chrom'
  gr[["score"]] # unlists over the chromosomes
  score(gr)
  gr[1][["score"]] # equal to score[1:3]
}
\keyword{classes}
\keyword{methods}