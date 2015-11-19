\name{CD.solver}
\alias{CD.solver}
\alias{CD.solver2}
\alias{CD.solver3}
\alias{CD.solver4}
\alias{ad2xyang}
\alias{single_image_resid}
\alias{coor_resid}
\alias{imageCD_solver}
\title{
Measure central world coordinates and rotation angles
}
\description{
Measures central world coordinates and rotation angles for each frame 
in a data set. Allows options to specify which stars to use, select a 
lower SNR limit, have multiple core processing, provide independent star
coordinate guesses, and provide independent central coordinates and 
rotation angle starting guesses.
}
\usage{
CD.solver(data, realz = NA, cores = 4, SNR_cut = 20, 
star_id = NA, cda = NA)

CD.solver2(data, realz = NA, cores = 4, SNR_cut = 20, 
star_id = NA, cda = NA)

CD.solver3(data, realz = NA, cores = 4, SNR_cut = 20, 
star_id = NA, cda = NA)

CD.solver4(data, realz = NA, cores = 4, SNR_cut = 20, 
star_id = NA, cda = NA)
}
\arguments{
  \item{data}{
Dataset containing all the individual detections
(i.e., Output of read.in.data).
}
  \item{realz}{
Optional initial guesses for world coordinates of stars. Should have same 
length as index$star_id
}
  \item{cores}{
Number of cores to run this function for multicore processing. 
Defaults to 4.
}
  \item{SNR_cut}{
Specify which detections should be ignored based on SNR lower limit. Defaults
to SNR < 20.
}
  \item{star_id}{
Option to specify which stars to use in calculations. All other detections
are ignored. Defaults to using all available stars.
}
  \item{cda}{
Optional initial guesses for individual image world coordinates and 
rotation. Should have same length as data$image_key$image_id
}
}
\details{
Since true coordinates of stars are not initially  known, the function relies
on an iterative approach.
Unless cda is specified, the first guess is derived from the image headers
recorded in data$image_key. 
Unless realz is specified, the estimate for star world coordinates comes 
from calc_all (We strongly recommend using 2MASS coordinates for realz).
Function then iterates measuring image parameters using a least squares
approach and recalculating star world coordinates. 
Outlier rejection is also applied twice. 
Uses internal functions: single_image_resid, coor_resid, and imageCD_solver,
which are not designed for external computations. The time-dependent scale
factor is estimated by linear interpolation between the known points
already measured. 

CD.solver2 does not revaluate the world coordinates but rather
matches detections to 
realz after applying aggressive outlier rejection multiple times. 
This option should be only used with an input IRAC catalog from a different
dataset that has been solved with CD.solver. Or the values of realz are 
very accurate. 

CD.solver3 follows the algorithm of CD.solver after first estimating 
the time-dependent scale factor. This scale factor should be measured 
for data observed at HMJDs greater than 3 months away from the measured
scale factors (see code for dis_scale for a list of known HMJDs).

CD.solver4 is CD.solver2 with an time-dependent scale factor estimation
routine first. 
}
\value{
Returns a matrix with columns defined as RA, DEC, angle, nums, 
and Residuals.

RA and DEC are both in degrees, angle in radians, and residuals in the 
average arcsecond residual between mean position of star and an individual 
detection. 'nums' is an integer showing the number of stars used
to measure the other three parameters. 

If CD.solver3 or CD.solver4 is used a list is returned with the 
first component being the previously described matrix. The second
component is a vector with the values of the estimated scale factor
for [3.6] and [4.5]. All other functions that take the output of CD.solver
are constructed to be able to take either a matrix or a list.
}
\author{
Taran Esplin
}
\examples{
\dontrun{
#This takes several minutes even with multiple cores so don't casually run
data(data1)
CDtest <- CD.solver(data1)
#compare above results with 
data(CD1)
#where a realz was used in CD.solver2
}
}
\keyword{ methods }
