\name{coor.calc}
\alias{coor.calc}
\alias{pix_correc}
\alias{dis_correc}
\alias{dis_scale}
\alias{apply_poly2d_simp}
\alias{dimvec2d}
\alias{ca_pix1}
\alias{wa_pars1}
\alias{wa_pars2}
\title{
Calculate World Coordinates from Pixel Coordinates
}
\description{
Basic calculator of world coordinates from pixel coordinates using a 
specified pixel bias corrector and 
optical distortion correction given the central
world coordinates of the chips central pixel and rotation angle. Generally
these latter variables will be measured by CD.solver.
}
\usage{
coor.calc(pix, pars, crsa, u, v, data, ch)
}
\arguments{
  \item{pix}{
Pixel bias correction polynomial - only cryogenic [3.6] 
observations need this parameter (ca_pix1 for [3.6]). 
If calculating post cryogenic data coordinates enter NA
}
  \item{pars}{
Optical distortion correction polynomial: [3.6] = wa_pars1; 
[4.5] = wa_pars2
}
  \item{crsa}{
Vector of world coordinates of the chips central pixel and rotation angle
c(RA, DEC, angle). If calculating world coordinates for multiple detections,
then rbind the respective c(RA, DEC, angle) for each detection.
}
  \item{u}{
x pixel coordinate with central pixel set to 0 (i.e., a range of -128 to
128)
}
  \item{v}{
y pixel coordinate with central pixel set to 0 (i.e., a range of -128 to
128)
}
  \item{data}{
output of read.in.data for the current dataset or a 2 element vector.
}
  \item{ch}{
integer value of which channel coordinates are being calculated for:
[3.6] = 1; [4.5] = 2.
}
}
\details{
Function first corrects for pixel bias, then optical distortion,
and finally projects the pixel coordinates into world coordinates using
the gnomonic projection. 
Other functions linked to this help page are internal 
to the coor.calc and should not be used by the user.

If data is the output of read.in.data, the function estimates
the time-dependent scale factor using the observation dates recorded
in the image headers. 
If data is a two element vector it assumes that the two numbers
are the [3.6] and [4.5] scale-factor-correcting ratios measured 
independently (from CD.solver3 or CD.solver4).  
}
\value{
Returns matrix of calculated RA and DEC (alpha and delta)
}
\author{
Taran Esplin
}
\examples{
data(CD1,ca_pix1,wa_pars1,data1)
options(digits=10)
#using a measured scale factor
coor.calc(ca_pix1,wa_pars1,CD1[[1]][1,],-100,104,CD1[[2]],1)
#estimating a scale factor from HMJD.
coor.calc(ca_pix1,wa_pars1,CD1[[1]][1,],-100,104,data1,1)
#the difference for this point in the array is ~2 mas
}
\keyword{ methods }
