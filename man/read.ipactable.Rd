\name{read.ipactable}
\alias{read.ipactable}
\title{
Read in Table output from IPAC
}
\description{
Tables from NASA/IPAC Infrared Science Archive have a specific format that
is not read into R through one of the normal table readers. This function 
will read an IPAC table and identify the columns.
}
\usage{
read.ipactable(file1)
}
\arguments{
  \item{file1}{
Character string containing pathway to file that needs to be read.
}
}
\value{
Returns a dataframe with the columns and column classes indicated in the
header to the IPAC table. 
}
\author{
Taran Esplin
}
\examples{
\dontrun{
	#This function cannot be easily provided with an example
	}
}
\keyword{ methods }
