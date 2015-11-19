\name{calc.all.ch12}
\alias{calc.all.ch12}
\alias{calc_all}
\title{
Calculate the 'average' RA and DEC for objects with multiple detections
}
\description{
Calculate 'average' RA and DEC of all sources from a full dataset of
detections from one epoch. calc_all is internal for CD.solver but can
be used for actual data, though not recommended.
}
\usage{
calc.all.ch12(data, cda, cores = 4, whichstars = NA, outlier = TRUE)
calc_all(pix1, pars1, pix2, pars2, cda, data, whichstars = NA, 
outlier = TRUE, cores = 4)
}
\arguments{
  \item{pix1}{
Pixel bias correction polynomial: Cryogenic = ca_pix1;
Post Cryogenic = NA
}
  \item{pars1}{
Optical distortion correction polynomial: [3.6] = wa_pars1
}
  \item{pix2}{
Pixel bias correction polynomial: Cryogenic = NA;
Post Cryogenic = NA
}
  \item{pars2}{
Optical distortion correction polynomial: [4.5] = wa_pars2
}
  \item{cda}{
Matrix containing c(RA,DEC,angle) for each frame in the dataset. 
Should be the same length as data$image_key$image_id. (The output 
of CD.solver or CD.solver2). Can also be a list containing the matrix
and a vector with an estimated scale factor (i.e., the output of
CD.solver3 or CD.solver4).
}
  \item{data}{
Dataset containing all the individual detections
(i.e., Output of read.in.data).
}
  \item{whichstars}{
Vector containing the star numbers that should be evaluated. If NA, 
program will calculate all stars coordinates. 
}
  \item{outlier}{
Set to FALSE if no outlier rejection desired.
This parameter has no effect calc.all.ch12.
}
  \item{cores}{
Number of cores to run this function for multicore processing. 
Defaults to 4.
}
}
\details{
calc_all calculates the straight mean of individual detections, 
optionally without outlier rejection (outlier=FALSE).
This version is used inside CD.solver but should generally not be used 
for full datasets. 

calc.all.ch12 computes the coordinates using a 
algorithm depending on the number of detections. 
For one detection, that detections coordinates are returned.
For 2 or more detections the median and mad (median absolute 
deviation) are returned. If an object
has more the 4 detections in [3.6] and the mad([3.6] detections)*2 < 
mad(all detections) then all [4.5] detections are ignored.
A similar test is applied to the [4.5] detections. 
}
\value{
calc.all.ch12 returns data.frame with following components:
\item{RA}{Calculated RA}
\item{DEC}{Calculated DEC}
\item{e_RA}{MAD of RA}
\item{e_DEC}{MAD of DEC}
\item{n_detections}{number of detections used to calculate coordinates}
\item{ch1_e_RA}{MAD of RA for just [3.6] detections}
\item{ch1_e_DEC}{MAD of DEC for just [3.6] detections}
\item{ch1_n_detections}{total number of detections from [3.6]}
\item{ch1_FLUX}{median FLUX for [3.6]}
\item{ch1_SNR}{median SNR for [3.6]}
\item{ch2_e_RA}{MAD of RA for just [4.5] detections}
\item{ch2_e_DEC}{MAD of DEC for just [4.5] detections}
\item{ch2_n_detections}{total number of detections from [4.5]}
\item{ch2_FLUX}{median FLUX for [4.5]}
\item{ch2_SNR}{median SNR for [4.5]}
calc_all returns matrix with columns of (mean RA, mean DEC, standard
deviation RA, standard deviation DEC, total number of detections). It does
not distinguish [3.6] and [4.5] detections. If outlier is set to TRUE,
a 2-sigma clip is made for sources with more than 2 detections and values
are then recomputed.
}
\author{
Taran Esplin
}
\examples{
data(data1,CD1,wa_pars1,wa_pars2,ca_pix1)
options(digits=10)
t1 <- calc.all.ch12(data1,CD1,1,c(1542,1840,1844,1847,2074,2089,2133))
t1 <- t1[!is.na(t1[,1]),]
t2 <- calc_all(ca_pix1,wa_pars1,NA,wa_pars2, 
  CD1,data1,c(1542,1840,1844,1847,2074,2089,2133),TRUE,1)
t2 <- t2[!is.na(t2[,1]),]
t1
t2
}
\keyword{ methods }
