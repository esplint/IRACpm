\name{mucalc}
\alias{mucalc}
\title{
Calculate relative proper motions of objects over multiple epochs
}
\description{
Measures the apparent motion of stars between multiple epochs using
a linear model of the motion. Optional use of the error on individual 
measurements normalized by the number of detections as weights. 
}
\usage{
mucalc(data, year, weight = FALSE)
}
\arguments{
  \item{data}{
Array containing all coordinates organized by epoch. 
}
  \item{year}{
Vector containing the time of each epoch in years. 
}
  \item{weight}{
Logical: if the errors on the coordinates divided by the number of detections 
should be used as weights in
the linear model fit. Defaults to FALSE.
}
}
\details{
The structure of data should be an array with 3 dimensions. The first 
dimension should be the named star_id and the length of star_id 
(see output of \code{\link{read.in.data}}).
The second dimension should be the epochs labeled ('1','2','3'...).
And the final dimension labeled 'RA', 'DEC', 'e_RA', 'e_DEC', 'n_detections'
(see \code{\link{calc.all.ch12}}). 

This function will not work if any coordinate has no associated error.
}
\value{
Returns matrix with columns labeled as mu.ra, mu.dec, mu.ra.sig,
and mu.dec.sig, which correspond to the apparent proper motion in the 
RA and DEC direction, and their respective errors. Row names are the star
id as indicated in the index object.
}

\author{
Taran Esplin
}

\examples{
#This takes a while to run.
data(epochs3)
str(epochs3)
test <- mucalc(epochs3,c(0.56,5.10,9.77))
test <- test[!is.na(test[,1]),]
test[1:10,]
}
\keyword{ methods }

