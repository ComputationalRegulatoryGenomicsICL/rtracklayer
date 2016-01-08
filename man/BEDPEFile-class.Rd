\name{BEDPEFile-class}
\docType{class}

%% Classes:
\alias{class:BEDPEFile}

%% Constructor:
\alias{BEDPEFile}

%% Import:
\alias{import,BEDPEFile,ANY,ANY-method}
\alias{import.bedpe}
\alias{import.bedpe,ANY-method}

%% Export:
\alias{export,ANY,BEDPEFile,ANY-method}
\alias{export,GInteractions,BEDPEFile,ANY-method}
\alias{export.bedpe}
\alias{export.bedpe,ANY-method}

\title{BEDPEFile objects}

\description{
  These functions support the import and export of the BEDPE format.
}

\usage{
\S4method{import}{BEDPEFile,ANY,ANY}(con, format, text,
                   genome = NA, seqinfo = NULL, extraCols = character())
import.bedpe(con, ...)

\S4method{export}{ANY,BEDPEFile,ANY}(object, con, format, ...)
\S4method{export}{GInteractions,BEDFile,ANY}(object, con, format,
                  append = FALSE, ignore.strand = FALSE)
export.bedpe(object, con, ...)
}

\arguments{
  \item{con}{A path, URL, connection or \code{BEDPEFile} object. For the
    functions ending in \code{.bedpe} the file format is indicated by the
    function name. For the base \code{export} and \code{import} functions, the
    format must be indicated another way. If \code{con} is a path, URL or
    connection, either the file extension or the \code{format} argument needs
    to be \dQuote{bedpe}. Compressed files (\dQuote{gz}, \dQuote{bz2} and
    \dQuote{xz}) are handled transparently.
  }
  \item{object}{The object to export, should be a \code{GInteractions} or
    something coercible to a \code{GInteractions}.
  }
  \item{format}{If not missing, should be \dQuote{bedpe}.}  
  \item{text}{If \code{con} is missing, a character vector to use as the input.}
  \item{genome}{The identifier of a genome, or \code{NA} if
    unknown. Typically, this is a UCSC identifier like \dQuote{hg19}. An
    attempt will be made to derive the \code{seqinfo} on the return
    value using either an installed BSgenome package or UCSC, if network
    access is available.
  }
  \item{ignore.strand}{Whether to output the strand when not required
    (by the existence of later fields).
  }
  \item{seqinfo}{If not \code{NULL}, the \code{Seqinfo} object to set on
    the result. If the \code{genome} argument is not \code{NA}, it must
    agree with \code{genome(seqinfo)}.
  }
  \item{extraCols}{A character vector in the same form as
    \code{colClasses} from \code{\link{read.table}}.  It should indicate
    the name and class of each extra/special column to read from the BEDPE
    file. As BEDPE does not encode column names, these are assumed to be
    the last columns in the file. If extraCols is not specified and extra 
    columns are present, they will be read but without names or type coercion.
  }
  \item{append}{If \code{TRUE}, and \code{con} points to a file path,
    the data is appended to the file. Obviously, if \code{con} is a
    connection, the data is always appended.
  }
  \item{...}{Arguments to pass down to methods to other methods.  }
}

\value{
  A \code{GInteractions} with the metadata columns described in the details.
}

\details{
  The BEDPE format is a tab-separated table of interval pairs, with annotations
  like name, score. It is not an official UCSC-supported format, and
  is incompatible with other BED and BED-like formats, but allows
  representation of paired intervals between chromosomes.

  To visualise interaction data in UCSC, export to BED12 is
  reccomended, which allows visualisation but can cause loss of data
  for trans-chromosomal and stranded interactions.

  The following is the mapping of BEDPE elements to a \code{GInteractions}
  object. NA values are allowed only where indicated. These appear as
  a \dQuote{.} in the file. The first ten columns are required.
  Extra columns can only be included to the right of these columns.
  Upon export, default values are used to automatically pad the table,
  if necessary.
  \describe{
    \item{chrom1, start1, end1, chrom2, start2, end2}{the \code{ranges} 
        components for the starts/ends of interactions. In the
        internal representation, these are compressed into one
        \code{GRanges} object in order to save memory.
    }
    \item{name}{character vector (NA's allowed) in the \code{name}
      column; defaults to NA on export. 
    }
    \item{score}{numeric vector in the \code{score}
      column, accessible via the \code{$} accessor. Defaults to 0
      on export.
    }
    \item{strand1, strand2}{strand factor (NA's allowed) for the \code{strand}
      column in the anchors; defaults to NA on export.
    }
    \item{Additional user-defined fields}{These are passed from the
        \code{mcols} of the \code{GInteractions} object. Be aware that
        column names are not imported/exported, however many parsers
        (rtracklayer::import included) will allow comments preceeded
        by \code{#} at the top of a file.
    }
  }
}

\section{BEDPEFile objects}{
  The \code{BEDPEFile} class extends \code{\linkS4class{RTLFile}} and is a
  formal represention of a resource in the BEDPE format.
  To cast a path, URL or connection to a \code{BEDPEFile}, pass it to
  the \code{BEDPEFile} constructor. 
}

\author{Malcolm Perry, Michael Lawrence}

\references{
  \url{http://bedtools.readthedocs.org/en/latest/content/general-usage.html}
}

\examples{
  test_path <- system.file("tests", package = "rtracklayer")
  test_bed <- file.path(test_path, "test.bedpe")

  test <- import(test_bedpe)
  test
  
  test_bedpe_file <- BEDPEFile(test_bedpe)
  import(test_bedpe_file)
  
  test_bedpe_con <- file(test_bedpe)
  import(test_bedpe_con, format = "bedpe")
  close(test_bedpe_con)
  
  import(test_bedpe, genome = "hg19")

  extraCols = c("numeric")
  namse(extraCols) = "FDR"
  import(test_bedpe, extraCols = extraCols)

\dontrun{
  test_bedpe_out <- file.path(tempdir(), "test.bedpe")
  export(test, test_bedpe_out)
  
  test_bedpe_out_file <- BEDPEFile(test_bedpe_out)
  export(test, test_bedpe_out_file)

  export(test, test_bedpe_out, name = "Alternative name")
  
  test_bedpe_gz <- paste(test_bedpe_out, ".gz", sep = "")
  export(test, test_bedpe_gz)
  
  # bedpe_text <- export(test, format = "bedpe")
  # test <- import(format = "bedpe", text = bedpe_text)
}
}

\keyword{methods}
\keyword{classes}
\keyword{InteractionSet}
\keyword{GInteractions}
