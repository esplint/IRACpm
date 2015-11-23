# IRACpm
IRACpm R Package: Applies a 7-8 order distortion correction to IRAC astrometric data from the Spitzer Space Telescope
and includes a function for measuring apparent proper motions between different Epochs.

#Instructions for installation into R

install.packages("devtools")

library(devtools)

install_github("esplint/IRACpm")

#Basic work flow:

1) Read in files containing output from the Spitzer Science Center’s APEX single frame module
form MOPEX using read.in.data.

2) Measure image central world coordinates and rotations with CD.solver, CD.solver2, CD.solver3,
or CD.solver4

3) Calculate average coordinates for each star of interest with calc.all.ch12

4) Repeat for other epochs

5) Run mucalc to measure apparent proper motions
(If accurate relative astrometry is wanted without proper motions, just follow steps 1-3.)

#Example

Example datasets and output for CD.solver is CD1, 
read.in.data is data1, and input data for mucalc is epochs3

To just convert pixel coordinates to World Coordinates using the distortion corrections measured
follow the example listed below.

data(CD1,ca_pix1,wa_pars1,data1)

options(digits=10)

--using a measured scale factor
coor.calc(ca_pix1,wa_pars1,CD1[[1]][1,],-100,104,CD1[[2]],1)
#estimating a scale factor from HMJD.
coor.calc(ca_pix1,wa_pars1,CD1[[1]][1,],-100,104,data1,1)
#the difference for this point in the array is ~2 mas
