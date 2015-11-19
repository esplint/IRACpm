\name{align}
\alias{align}
\alias{correct_set}
\alias{select_bright}
\alias{cut_error}
\alias{E1}
\alias{E2}
\title{
Align One Epoch with Another
}
\description{
Aligns output of calc_all_ch12 of one epoch to another with a simple 
median ra and dec shift. This can only correct small offset.  
}
\usage{
align(data1, data2, ra, dec, exclude, ERROR=75, SNR=30)
}

\arguments{
  \item{data1}{
Output of calc_all_ch12 for the epoch to be shifted.
}
  \item{data2}{
Output of calc_all_ch12 for the epoch that is the reference.
}
  \item{ra}{
	two component vector for RA boundaries. Will only shift sources
within this bound
  }
  \item{dec}{
  	two component vector for DEC boundaries. Will only shift sources
within this bound
  }
  \item{exclude}{
  	logical with length(data1[,1]). TRUE for sources that have known
  	large motions between epochs.
  }
  \item{ERROR}{
  	Sources with errors below this limit are also exclude from calculations 
  	of offsets. Default 75 mas.
  }
  \item{SNR}{
  	Sources with SNR below this limit are also exclude from calculations 
  	of offsets. Default 30.
  }
}
\details{
Shifts with the following steps:
(1) Identifies sources that are not excluded manually or with ERROR/SNR.
(2) Identifies all sources that fit the RA and DEC boundaries.
(3) Makes a two sigma clip for large motions 
(4) shifts data1 to align with data2 with independent RA and DEC shifts.
Uses the internal functions correct_set, select_bright, and cut_error.
}
\value{
Returns data1 with RA and DEC shifted.
}
\author{
Taran Esplin
}
\examples{
data(E1,E2)
#plot((E2[,1:2]-E1[,1:2])*3600,pch=20,cex=.1,xlim=c(-0.3,0.3),
#  ylim=c(-0.3,0.3))
E2n <- align(E2,E1,c(48,55),c(21,35),rep(FALSE,length(E1[,1])),75,15)
#dev.new()
#plot((E2n[,1:2]-E2[,1:2])*3600,pch=20,cex=.1,xlim=c(-0.05,0.05),
#  ylim=c(-0.05,0.05))
}
\keyword{ methods }
