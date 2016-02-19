library('foreach')
library('doMC')
library('astro')
library('R.utils')

globalVariables(c('wa_pars1','wa_pars2','ca_pix1','j','badpixch1','badpixch2'))

read.in.data <- function(data.ch1=NA,data.ch2=NA,ch1.imag=NA,ch2.imag=NA,star_id,cores=4) {
#makes and fills a table with all the data from the orginal files:
#after running this once you don't need to read the files again.
#data to read in should be a vector of full path names and images full path.

#in the individual files should have the following format:
#'RA','DEC','id','ra.calc','dec.calc','x','y','mag','snr','sep'
#generally made by topcat matching your files to an allstar list.

  if (!is.na(data.ch1[1])) {
    #count the number of lines of all files and make a huge matrix.
    ch1.length <- length(data.ch1)
    ch1.lens <- rep(NA,ch1.length)
    for (j in 1:ch1.length) ch1.lens[j] <- countLines(data.ch1[j])
  } else {
    ch1.length <- 0
    ch1.lens <- 0
  }
  if (!is.na(data.ch2[1])) {
    ch2.length <- length(data.ch2)
    ch2.lens <- rep(NA,ch2.length)
    for (j in 1:ch2.length) ch2.lens[j] <- countLines(data.ch2[j])
  } else {
    ch2.length <- 0
    ch2.lens <- 0
  }


  #data frame for individual data points
  gen_set <- data.frame(matrix(nrow=sum(ch2.lens)+sum(ch1.lens),ncol=6, dimnames=list(NULL,c('image_id','star_id','x','y','FLUX','SNR'))))

  #data for each image in the data set
  imkey <- data.frame(matrix(nrow=ch2.length+ch1.length,ncol=10,dimnames=list(NULL,c('image_id','CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2','AORKEY','HMJD','EXPTIME'))))

  oldspot <- 0 #counter of place in gen_set

#fill in Channel 1 frames:
  if (!is.na(data.ch1[1])) {
    for (j in 1:length(data.ch1)) {
      temp <- read.table(data.ch1[j],skip=1,col.names=c('RA','DEC','id','ra.calc','dec.calc','x','y','mag','snr','sep'),na.strings='""')
      #shift to reference from the central pixels
      temp$x <- temp$x - 128 ; temp$y <- temp$y - 128

      #identifying stars for 'sid' which is star id
      ss2 <- as.numeric(temp$id)
 
      if (!is.na(ch1.imag[1])) {
        imt <- read.fitshdr(ch1.imag[j])
        exp <- as.numeric(imt[which(imt[,1] == "EXPTIME"),2])
        hjd <- as.numeric(imt[which(imt[,1] == "HMJD_OBS"),2])
        aor <- as.numeric(imt[which(imt[,1] == "AORKEY"),2])
        cra <- as.numeric(imt[which(imt[,1] == "CRVAL1"),2])
        cdec <- as.numeric(imt[which(imt[,1] == "CRVAL2"),2])
        cd11 <- as.numeric(imt[which(imt[,1] == "CD1_1"),2]) 
        cd12 <- as.numeric(imt[which(imt[,1] == "CD1_2"),2])
        cd21 <- as.numeric(imt[which(imt[,1] == "CD2_1"),2])
        cd22 <- as.numeric(imt[which(imt[,1] == "CD2_2"),2])
      } else {
        exp <- NA
        hjd <- NA
        aor <- NA
        cra <- NA
        cdec <- NA
        cd11 <- NA
        cd12 <- NA
        cd21 <- NA
        cd22 <- NA
      } 

      newspot <- oldspot + length(ss2)
      gen_set[(oldspot+1):newspot,] <- cbind(rep(j,length(ss2)),ss2,temp$x,temp$y,temp$mag,temp$snr)
      

      imkey[j,] <- cbind(j,cra,cdec,cd11,cd12,cd21,cd22,aor,hjd,exp)

      oldspot <- newspot
    }
  }

#Channel 2:
  if (!is.na(data.ch2[1])) {
    for (j in 1:length(data.ch2)) {
      temp <- read.table(data.ch2[j],skip=1,col.names=c('RA','DEC','id','ra.calc','dec.calc','x','y','mag','snr','sep'),na.strings='""')
      #shift to reference from the central pixels
      temp$x <- temp$x - 128 ; temp$y <- temp$y - 128

      #identifying stars for 'sid' which is star id
      ss2 <- as.numeric(temp$id)
  
      if (!is.na(ch2.imag[1])) {
        imt <- read.fitshdr(ch2.imag[j])
        exp <- as.numeric(imt[which(imt[,1] == "EXPTIME"),2])
        hjd <- as.numeric(imt[which(imt[,1] == "HMJD_OBS"),2])
        aor <- as.numeric(imt[which(imt[,1] == "AORKEY"),2])
        cra <- as.numeric(imt[which(imt[,1] == "CRVAL1"),2])
        cdec <- as.numeric(imt[which(imt[,1] == "CRVAL2"),2])
        cd11 <- as.numeric(imt[which(imt[,1] == "CD1_1"),2]) 
        cd12 <- as.numeric(imt[which(imt[,1] == "CD1_2"),2])
        cd21 <- as.numeric(imt[which(imt[,1] == "CD2_1"),2])
        cd22 <- as.numeric(imt[which(imt[,1] == "CD2_2"),2])
      } else {
        exp <- NA
        hjd <- NA
        aor <- NA
        cra <- NA
        cdec <- NA
        cd11 <- NA
        cd12 <- NA
        cd21 <- NA
        cd22 <- NA
      } 


      newspot <- oldspot + length(ss2)
      gen_set[(oldspot+1):newspot,] <- cbind(rep(-j,length(ss2)),ss2,temp$x,temp$y,temp$mag,temp$snr)
      imkey[j+ch1.length,] <- cbind(-j,cra,cdec,cd11,cd12,cd21,cd22,aor,hjd,exp)

      oldspot <- newspot
    }
  }

t1 <- !is.na(gen_set[,3])

imkey$AORKEY <- as.factor(imkey$AORKEY)

data <- list(data=gen_set[t1,],image_key=imkey)
data <- bad_pix_trim(data)
index <- star_index_maker(data,star_id,cores)

return(c(data,index))

}

bad_pix_trim <- function(data,edget=TRUE) {
  data(badpixch1,badpixch2)

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)

  bad <- numeric()
  for(j in 1:length(badpixch1$x1)) {
    
    bad <- c(bad,intersect(which(data$data$x > badpixch1$x1[j] & data$data$x < badpixch1$x2[j] & data$data$y > badpixch1$y1[j] & data$data$y < badpixch1$y2[j]),i1))
    
  }
  for(j in 1:length(badpixch2$x1)) {
    
    bad <- c(bad,intersect(which(data$data$x > badpixch2$x1[j] & data$data$x < badpixch2$x2[j] & data$data$y > badpixch2$y1[j] & data$data$y < badpixch2$y2[j]),i2))
    
  }

    edge1 <- intersect(which(data$data$x < -115 | data$data$x > 123 | data$data$y < -122 | data$data$y > 125),i1)
    
    edge2 <- intersect(which(data$data$x < -120 | data$data$x > 124 | data$data$y < -122 | data$data$y > 125),i2)

    abd <- c(bad,edge1,edge2)


  data$data$x[abd] <- NA
  data$data$y[abd] <- NA

  #eliminate any sources with SNR < 0 
  bbad <- which(data$data$SNR < 0)
  data$data$x[bbad] <- NA
  data$data$y[bbad] <- NA
  data$data$star_id[bbad] <- NA
  data$data$SNR[bbad] <- NA


  #eliminate any sources above our measured flux limit.... 
  #flux limits: 728000 FLUX * SEC exposure CH1
  #821600 * SEC exposure for CH2
  badf <- foreach (j=1:length(data$image_key$image_id),.combine=c) %dopar% {
    if (data$image_key$image_id[j] > 0) {
      which(data$data$FLUX > 728000 / data$image_key$EXPTIME[j] & data$data$image_id == data$image_key$image_id[j])
    } else {
      which(data$data$FLUX > 821600 / data$image_key$EXPTIME[j] & data$data$image_id == data$image_key$image_id[j])
    }
  }
  
  if (length(badf) != 0) {data$data[badf,3:4] <- NA}

return(data)
}

star_index_maker <- function(data,ids,cores=4) {
#this will give the indices of all detections of each allstars star...
#get rid of single detections....
  stars <- length(ids)
  registerDoMC(cores)

  
  if (length(stars) > 50000) {
    index <- list()
    start <- 1
    stop <- 50000
    for (ss in 1:ceiling(stars/50000)) {
      index.nex <- foreach (j=start:stop) %dopar% {
        t1 <- which(data[[1]][,2] == ids[j])
        if (length(t1) > 0) { 
          t1
        } else {
          NA
        }
      }
      index <- c(index,index.nex)
      start <- stop + 1
      stop <- min(c(stop+50000,stars))
    }
  } else {
    index <- foreach (j=1:stars) %dopar% {
      t1 <- which(data[[1]][,2] == ids[j])
      if (length(t1) > 0) { 
        t1
      } else {
        NA
      }    
    }
  }
  registerDoMC(1)


  nde <- rep(NA,length(index))
  for (j in 1:length(index)) {
    t1 <- length(index[[j]])
    if (is.na(index[[j]][1])) {
      nde[j] <- 0 
    } else {
      nde[j] <- t1
    }
  } 

return(list(index=index,n_detections=nde,star_id=ids))
}

read.ipactable <- function(file1) {

  con <- file(file1)
  open(con)
  spot <- 0
  while (substr(line <- readLines(con, n=1, warn=F),1,1) == '\\') {
    spot <- spot + 1
  }
  tdh <- line

  hdr <- gsub("^\\s+|\\s+$", "", strsplit(tdh,'|',fixed=T)[[1]])
  hdr <- hdr[2:length(hdr)]

  line <- readLines(con, n=1, warn=F) 
  units <- gsub("^\\s+|\\s+$", "", strsplit(line,'|',fixed=T)[[1]])
  units <- units[2:length(units)]

  units[units == 'long'] <- 'integer'
  units[units == 'int'] <- 'integer' 
  units[units == 'char'] <- 'character'
  units[units == 'double'] <- 'numeric'
  units[units == 'float'] <- 'numeric'

  units[units == 'i'] <- 'integer'
  units[units == 'c'] <- 'character'
  units[units == 'd'] <- 'numeric'
  units[units == 'r'] <- 'numeric'


  spot <- spot + 4
  close(con)

  tab <- read.table(file1,skip=spot,col.names=hdr,colClasses=units,na.strings=c('null','-'))
  return(tab)
}


CD.solver <- function(data,realz=NA,cores=4,SNR_cut=20,star_id=NA,cda=NA) {
data(wa_pars1,wa_pars2,ca_pix1)
  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  c1 <- which(data$image_key$image_id > 0)
  c2 <- which(data$image_key$image_id < 0)

  if (mean(data$image_key$HMJD,na.rm=T) > 54976) {WARM <- TRUE} else {WARM <- FALSE}

  if (length(cda) == 2) {
    sdata <- cda[[2]]
    cda <- cda[[1]]
  } else {
    sdata <- data
  }

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,sdata,1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,sdata,2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,sdata,1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,sdata,2)   
    }
  }
  
  if (is.na(star_id[1])) {
    star_id <- data$star_id
  }
#number of stars:
	nst <- length(data$index)
#number of frames
	nfr <- length(data$image_key$image_id)
  nro <- length(data[[1]]$image_id)
#median SNRs
  nde2 <- which(data$n_detections != 0)
  snr <- rep(NA,length(data$index))
  for (j in 1:length(nde2)) {
    snr[nde2[j]] <- median(data$data$SNR[data$index[[nde2[j]]]],na.rm=T)
  }

	
  nde2 <- which(data$n_detections > 2 & snr > SNR_cut)
  nde2 <- intersect(nde2,match(star_id,data$star_id))
  star_id <- data$star_id


  if (is.na(cda[1])) {
        cda <- matrix(0, nrow = nfr, ncol = 3)
        cda[, 1:2] <- as.matrix(data[[2]][, 2:3])
        cda[, 3] <- atan2(-data[[2]][,5],data[[2]][,7])-pi
    }	


#intial guess at coordinates:
	if (length(realz) == 1) {
		i.coor1 <- calc_all(pix1,pars1,pix2,pars2,cda,data,data$star_id[nde2],FALSE,cores) 
	} else {
		i.coor1 <- matrix(NA,nrow=length(data$star_id),ncol=2)
    i.coor1[nde2,] <- as.matrix(realz[nde2,1:2])
	}

	cda1 <- imageCD_solver(data,i.coor1,cda,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda1) == 4) { cda1 <- t(cda1)}
  nums <- cda1[,4]
  cda1 <- cda1[,1:3]
  if (length(cda1) == 3) { cda1 <- t(cda1)}

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda1[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda1[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2


  rre <- i.coor1[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])
  s2 <- RMS(rresd[,2])

	k1 <- which(abs(rresd[,1]) > 6*s1 | abs(rresd[,2]) > 6*s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor2 <- calc_all(pix1,pars1,pix2,pars2,cda1,data,data$star_id[nde2],TRUE,cores)

	cda2 <- imageCD_solver(data,i.coor2,cda1,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}
  #res2 <- single_image_resid(data,i.coor2,cda2,pix1,pars1,pix2,pars2,cores,star_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda2[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda2[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- i.coor2[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])
  s2 <- RMS(rresd[,2])

	k1 <- which(abs(rresd[,1]) > 6*s1 | abs(rresd[,2]) > 6*s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor3 <- calc_all(pix1,pars1,pix2,pars2,cda2,data,data$star_id[nde2],TRUE,cores)
	cda3 <- imageCD_solver(data,i.coor3,cda2,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda3) == 4) { cda3 <- t(cda3)}
  nums <- cda3[,4]
  cda3 <- cda3[,1:3]
  if (length(cda3) == 3) { cda3 <- t(cda3)}
  #res3 <- single_image_resid(data,i.coor3,cda3,pix1,pars1,pix2,pars2,cores,star_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda3[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda3[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- i.coor3[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])
  s2 <- RMS(rresd[,2])

	k1 <- which(abs(rresd[,1]) > 3*s1 | abs(rresd[,2]) > 3*s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor4 <- calc_all(pix1,pars1,pix2,pars2,cda3,data,data$star_id[nde2],TRUE,cores)
	cda4 <- imageCD_solver(data,i.coor4,cda3,pix1,pars1,pix2,pars2,cores,star_id)
  if (length(cda4) == 4) { cda4 <- t(cda4)}
  nums <- cda4[,4]
  cda4 <- cda4[,1:3]
  if (length(cda4) == 3) { cda4 <- t(cda4)}
  res4 <- single_image_resid(data,i.coor4,cda4,pix1,pars1,pix2,pars2,cores,star_id)

	return(cbind(cda4,nums,res4))
}

calc_all <- function(pix1,pars1,pix2,pars2,cda,data,whichstars=NA,outlier=TRUE,cores=4) {
  #which are ch1 and ch2

  if (is.na(whichstars[1])) { 
    nde <- which(data$n_detections > 0)
  } else {
    nde <- match(whichstars,data$star_id)
  }

  if (length(cda) == 2) {
    sdata <- cda[[2]]
    cda <- cda[[1]]
  } else {
    sdata <- data
  }

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  nro <- length(data[[1]]$image_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)  
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  mens <- matrix(NA,nrow=length(data$index),ncol=5)
  registerDoMC(cores) 
  mens1 <- foreach (j=1:length(nde),.combine=rbind) %dopar% {
        if (length(data$index[[nde[j]]]) < 2) {
          if (!is.na(data$index[[nde[j]]])) {
            c(i.coor[data$index[[nde[j]]],],NA,NA,1)
          }
        } else {
          if ((ln <- length(data$index[[nde[j]]])) < 3 | !outlier) {
            c(colMeans(i.coor[data$index[[nde[j]]],],na.rm=T),apply(i.coor[data$index[[nde[j]]],],2,sd,na.rm=T),ln)
          } else {
            t1 <- i.coor[data$index[[nde[j]]],]
            tmen <- colMeans(i.coor[data$index[[nde[j]]],],na.rm=T)
            t1[,1] <- (t1[,1]-tmen[1])*cos(tmen[2]*pi/180)
            t1[,2] <- t1[,2]-tmen[2]
            sds <- apply(t1,2,sd,na.rm=T)
            gd <- which(t1[,1]^2 + t1[,2]^2 < (2*(sds[1]^2 + sds[2]^2)))
            c(colMeans(i.coor[data$index[[nde[j]]][gd],],na.rm=T),apply(i.coor[data$index[[nde[j]]][gd],],2,sd,na.rm=T),length(gd)  )
          }
        }
      }
  registerDoMC(1)
  mens[nde,] <- mens1   

return(mens)
}

coor.calc <- function(pix,pars,crsa,u,v,data,ch){
#gives the calculated values using my SIP compatable distoriton correction

  if (!is.na(pix[1])) {
    cor <- pix_correc(u,v,pix) 
    u <- cor[,1]
    v <- cor[,2]
  }

  if (!is.na(pars[1])) {
    discor <- dis_correc(u,v,pars,data,ch) 
    u <- discor[,1]
    v <- discor[,2]
  }

  if (length(crsa) == 3 | length(crsa) == 4 | length(crsa) == 5) {crsa <- t(crsa)}

  cra <- crsa[,1]*pi/180.
  cdec <- crsa[,2]*pi/180.
  xsi <- (cos(crsa[,3])*u - sin(crsa[,3])*v)*pi/180.
  eta <- (sin(crsa[,3])*u + cos(crsa[,3])*v)*pi/180.  
  phi <- atan2(xsi,-eta)
  theta <- atan(1/(sqrt(xsi^2 + eta^2)))
  alpha <- cra + atan2(-cos(theta)*sin(phi-pi),sin(theta)*cos(cdec) - cos(theta)*sin(cdec)*cos(phi-pi))
  delta <- asin(sin(theta)*sin(cdec) + cos(theta)*cos(cdec)*cos(phi-pi))
  return(cbind(alpha,delta)*180/pi)
}

pix_correc <- function(u,v,opti) {
  degx <- opti[[3]]
  degy <- opti[[4]]

  valx <- length(opti[[1]])
  valy <- length(opti[[2]])
  diff  <- 3:(30+1)
  diff2 <- diff
  for (j in 1:length(diff)) {diff2[j] <- sum(diff[1:j])}
  t1 <- which(valx - diff2 == 0)
  linearx <- F
  unitx <- F
  if (length(t1) == 0) { 
    t1 <- which(valx - c(diff2+2) == 0)
    linearx <- T
  }
  if (length(t1) == 0) { 
    t1 <- which(valx - c(diff2+3) == 0)
    linearx <- T
    unitx <- T
  }
  t1 <- which(valy - diff2 == 0)
  lineary <- F
  unity <- F
  if (length(t1) == 0) { 
    t1 <- which(valy - c(diff2+2) == 0)
    lineary <- T
  }
  if (length(t1) == 0) { 
    t1 <- which(valy - c(diff2+3) == 0)
    lineary <- T
    unity <- T
  }

  up <- u - floor(u - opti$cuts[1])
  vp <- v - floor(v - opti$cuts[2])

  keep <- which(!is.na(u) & !is.na(v))
  uc <- rep(NA,length(u))
  vc <- rep(NA,length(v))

  uc[keep] <- apply_poly2d_simp(up[keep],vp[keep],opti[[1]],degx,linearx,unitx)
  vc[keep] <- apply_poly2d_simp(up[keep],vp[keep],opti[[2]],degy,lineary,unity)


    uuc <- uc + u
    vvc <- vc + v

  return(cbind(uuc,vvc))

}

apply_poly2d_simp <- function(x,y,coef,deg,linear,unit) {
  ids <- dimvec2d(deg,linear,unit)
  cc <- rep(0,length(x))
  for(j in 1:length(coef)) {
    cc <- cc + coef[j] *x^ids[[1]][j] *y^ids[[2]][j]
  }
  return(cc)
}


dimvec2d <- function(deg,linear,unit) {
  xs <- rep(NA,length(sum(3:deg+1)))
  ys <- rep(NA,length(xs))
  for (dd in 2:deg) {
    xs[(sum(3:(dd+1)) - dd):sum(3:(dd+1))] <- dd:0
    ys[(sum(3:(dd+1)) - dd):sum(3:(dd+1))] <- 0:dd
  }
  if (linear) {
    xs <- c(1,0,xs)
    ys <- c(0,1,ys)
  }
  if (deg == 1) {
    xs <- c(1,0)
    ys <- c(0,1)
  }
  if (unit) {
    xs <- c(0,xs)
    ys <- c(0,ys)
  }
  if (deg == 0) {
    xs <- 0
    ys <- 0
  }
  return(list(xs,ys))
}


dis_scale <- function(data) {
  me <- mean(data$image_key$HMJD,na.rm=T)
  if (me > 56611) {
    sca1 <- 1.0000400
    sca2 <- 1.0000830
  }
  if (me > 56253 & me <= 56611) {
    sca1 <- (1.0000400-0.9998759)/(56611-56253) * (me - 56253) + 0.9998759
    sca2 <- (1.0000830-0.9998889)/(56611-56253) * (me - 56253) + 0.9998889
  }
  if (me > 55774 & me <= 56253) {
    sca1 <- (0.9998759-1)/(56253-55774) * (me - 56253) + 0.9998759
    sca2 <- (0.9998889-1)/(56253-55774) * (me - 56253) + 0.9998889
  }
  if (me > 55423 & me <= 55774) {
    sca1 <- (1-0.9999786)/(55774-55423) * (me - 55774) + 1
    sca2 <- (1-1.0000490)/(55774-55423) * (me - 55774) + 1
  }
  if (me > 55105 & me <= 55423) {
    sca1 <- (0.9999786-0.9999226)/(55423-55105) * (me - 55105) + 0.9999226
    sca2 <- (1.0000490-0.9999565)/(55423-55105) * (me - 55105) + 0.9999565
  }
  if (me > 54976 & me <= 55105) {
    sca1 <- 0.9999226
    sca2 <- 0.9999565
  }
  if (me > 54571 & me <= 54976) {
    sca1 <- 1.0000890
    sca2 <- 1.0003440
  }
  if (me > 53840 & me <= 54571) {
    sca1 <- (1.0000890 - 1.0001880)/(54571-53840) * (me - 53840) + 1.0001880
    sca2 <- (1.0003440 - 1.0003880)/(54571-53840) * (me - 53840) + 1.0003880
  }
  if (me > 53628 & me <= 53840) {
    sca1 <- (1.0001880-1.0000220)/(53840-53628) * (me - 53840) + 1.0001880
    sca2 <- (1.0003880-1.0002100)/(53840-53628) * (me - 53840) + 1.0003880
  }
  if (me > 53354 & me <= 53628) {
    sca1 <- (1.0000220-1.0000700)/(53628-53354) * (me - 53354) + 1.0000700
    sca2 <- (1.0002100-1.0002800)/(53628-53354) * (me - 53354) + 1.0002800
  }
  if (me > 52988 & me <= 53354) {
    sca1 <- (1.0000700-1.0000440)/(53354-52988) * (me - 53354) + 1.0000700
    sca2 <- (1.0002800-1.0003200)/(53354-52988) * (me - 53354) + 1.0002800
  }
  if (me <= 52988) {
    sca1 <- 1.0000440
    sca2 <- 1.0003200
  }
  return(c(sca1,sca2))
}


dis_correc <- function(u,v,opti,data,ch) {
  degx <- opti[[3]]
  degy <- opti[[4]]

  if (length(data) > 2) {
    sca <- dis_scale(data)
  } else {
    sca <- data
  }
  if (ch == 1) {
    opti[[1]] <- opti[[1]]*sca[1]
    opti[[2]] <- opti[[2]]*sca[1]
  } else {
    opti[[1]] <- opti[[1]]*sca[2]
    opti[[2]] <- opti[[2]]*sca[2]
  }

  valx <- length(opti[[1]])
  valy <- length(opti[[2]])
  diff  <- 3:(30+1)
  diff2 <- diff
  for (j in 1:length(diff)) {diff2[j] <- sum(diff[1:j])}
  t1 <- which(valx - diff2 == 0)
  linearx <- F
  unitx <- F
  if (length(t1) == 0) { 
    t1 <- which(valx - c(diff2+2) == 0)
    linearx <- T
  }
  if (length(t1) == 0) { 
    t1 <- which(valx - c(diff2+3) == 0)
    linearx <- T
    unitx <- T
  }
  t1 <- which(valy - diff2 == 0)
  lineary <- F
  unity <- F
  if (length(t1) == 0) { 
    t1 <- which(valy - c(diff2+2) == 0)
    lineary <- T
  }
  if (length(t1) == 0) { 
    t1 <- which(valy - c(diff2+3) == 0)
    lineary <- T
    unity <- T
  }

  keep <- which(!is.na(u) & !is.na(v))
  uc <- rep(NA,length(u))
  vc <- rep(NA,length(v))

 

  uc[keep] <- apply_poly2d_simp(u[keep],v[keep],opti[[1]],degx,linearx,unitx)
  vc[keep] <- apply_poly2d_simp(u[keep],v[keep],opti[[2]],degy,lineary,unity)


    uu <- uc
    vv <- vc

  return(cbind(uu,vv))
}

RMS <- function(num) sqrt(sum((num)^2,na.rm=T)/sum(!is.na(num)))

imageCD_solver <- function(data,realz,cd1,pix1,pars1,pix2,pars2,cores,star_id,outlier=T) {
#steps through each image and resolves for 
#central coordinates, with a known rotation angle...
#previous central coordinates are necessary:

  if (length(cd1) == 2) {
    sdata <- cd1[[2]]
    cda <- cd1[[1]]
  } else {
    sdata <- data
  }

  #number of images:
  inum <- length(data[[2]][,1])

  registerDoMC(cores)
  cd2 <- foreach (j=1:inum,.combine=rbind) %dopar% {
    t1 <- which(data[[1]]$image_id == data[[2]]$image_id[j])
    if (data[[2]]$image_id[j] > 0) { 
      crsa <- optim(cd1[j,1:3],fn=coor_resid,u=data[[1]]$x[t1],v=data[[1]]$y[t1],method='BFGS',rre=realz[match(data[[1]]$star_id[t1],star_id),1:2],pars=pars1,pix=pix1,data=sdata,ch=1)$par
      if (outlier) {
        calc.coords <- coor.calc(pix1,pars1,crsa,data[[1]]$x[t1],data[[1]]$y[t1],sdata,1)
        rre <- realz[match(data[[1]]$star_id[t1],star_id),1:2]
        res <- (rre[,1:2] - calc.coords)*cbind(cos(rre[,2]*pi/180),rep(1,length(calc.coords[,1])))*3600
        bd <- which(abs(res[,1]) > RMS(res[,1])*3 | abs(res[,2]) > RMS(res[,2])*3)
        t1 <- setdiff(t1,t1[bd])
        crsa <- optim(cd1[j,1:3],fn=coor_resid,u=data[[1]]$x[t1],v=data[[1]]$y[t1],method='BFGS',rre=realz[match(data[[1]]$star_id[t1],star_id),1:2],pars=pars1,pix=pix1,data=sdata,ch=1)$par
      }
    } else {
      crsa <- optim(cd1[j,1:3],fn=coor_resid,u=data[[1]]$x[t1],v=data[[1]]$y[t1],method='BFGS',rre=realz[match(data[[1]]$star_id[t1],star_id),1:2],pars=pars2,pix=pix2,data=sdata,ch=2)$par
      if (outlier) {
        calc.coords <- coor.calc(pix2,pars2,crsa,data[[1]]$x[t1],data[[1]]$y[t1],sdata,2)
        rre <- realz[match(data[[1]]$star_id[t1],star_id),1:2]
        res <- (rre[,1:2] - calc.coords)*cbind(cos(rre[,2]*pi/180),rep(1,length(calc.coords[,1])))*3600
        bd <- which(abs(res[,1]) > RMS(res[,1])*3 | abs(res[,2]) > RMS(res[,2])*3)
        t1 <- setdiff(t1,t1[bd])
        crsa <- optim(cd1[j,1:3],fn=coor_resid,u=data[[1]]$x[t1],v=data[[1]]$y[t1],method='BFGS',rre=realz[match(data[[1]]$star_id[t1],star_id),1:2],pars=pars2,pix=pix2,data=sdata,ch=2)$par
      }

    }
    c(crsa,sum(!is.na(realz[match(data[[1]]$star_id[t1],star_id),1])))
  }
  registerDoMC(1)
  return(cd2)
}

coor_resid <- function(crsa,pars,pix,u,v,rre,data,ch) {
  calc.coords <- coor.calc(pix,pars,crsa,u,v,data,ch)
  return(sum(((rre[,1:2] - calc.coords)*cbind(cos(rre[,2]*pi/180),rep(1,length(calc.coords[,1])))*3600)^2,na.rm=T) / sum(!is.na(rre[,1])))
}

single_image_resid <- function(data,realz,cds,pix1,pars1,pix2,pars2,cores,star_id) {
#individual image residuals from fit...
#for comparision:
  #resid <- rep(NA,length(cds[,1]))

  if (length(cds) == 2) {
    sdata <- cds[[2]]
    cda <- cds[[1]]
  } else {
    sdata <- data
  }

registerDoMC(cores)
#step over each image in the list and add the totals together...
  resid <- foreach (j=1:length(data[[2]][,1]),.combine=c) %dopar% {
    t1 <- which(data[[1]]$image_id == data[[2]][j,1])
    if (data[[2]][j,1] > 0) { 
      t2 <- coor_resid(as.numeric(cds[j,]),pars1,pix1,u=data[[1]]$x[t1],v=data[[1]]$y[t1],realz[match(data[[1]]$star_id[t1],star_id),1:2],sdata,1)
    } else {
      t2 <- coor_resid(as.numeric(cds[j,]),pars2,pix2,u=data[[1]]$x[t1],v=data[[1]]$y[t1],realz[match(data[[1]]$star_id[t1],star_id),1:2],sdata,2)
    }
    t2
  }
registerDoMC(1)
  return(resid)
}

calc.all.ch12 <- function(data,cda,cores=4,whichstars=NA,outlier=TRUE) {
data(wa_pars1,wa_pars2,ca_pix1)

  if (mean(data$image_key$HMJD,na.rm=T) > 54976) {WARM <- TRUE} else {WARM <- FALSE}

  if (length(cda) == 2) {
    sdata <- cda[[2]]
    cda <- cda[[1]]
  } else {
    sdata <- data
  }


  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- wa_pars1
    pars2 <- wa_pars2
  } else {
    pix1 <- ca_pix1
    pix2 <- NA
    pars1 <- wa_pars1
    pars2 <- wa_pars2
  }

  if (is.na(whichstars[1])) { 
    nde <- which(data$n_detections > 0)
  } else {
    nde <- rep(NA,length(whichstars))
    for (j in 1:length(nde)) {
      nde[j] <- which(data$star_id == whichstars[j])
    }
  }

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  nro <- length(data[[1]]$image_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  mens <- data.frame(matrix(NA,nrow=length(data$index),ncol=15,dimnames=list(NULL,c('RA','DEC','e_RA','e_DEC','n_detections','ch1_e_RA','ch1_e_DEC','ch1_n_detections','ch1_FLUX','ch1_SNR', 'ch2_e_RA','ch2_e_DEC','ch2_n_detections','ch2_FLUX','ch2_SNR'))))
  registerDoMC(cores) 
  mens1 <- foreach (j=1:length(nde),.combine=rbind) %dopar% {
      ii <- data$index[[nde[j]]]
      ii <- ii[which(!is.na(i.coor[ii,1]))]
      ln <- length(ii)
      if (ln < 2) {
        if (ln > 0) {
           c(i.coor[ii,],NA,NA,1,rep(NA,10))
        } else {
           c(NA,NA,NA,NA,0,rep(NA,10))
        }
      } else {
           if (ln < 4) {
            ln2 <- ln
            cor <- apply(i.coor[ii,],2,mean,na.rm=T)
            sd2 <- apply(i.coor[ii,],2,mad,na.rm=T)*c(cos(cor[2]*pi/180),1)
            tp <- intersect(ii,i1)
            bt <- intersect(ii,i2)           
            era1 <- mad(i.coor[tp,1],center=cor[1],na.rm=T)*cos(cor[2]*pi/180)
            era2 <- mad(i.coor[bt,1],center=cor[1],na.rm=T)*cos(cor[2]*pi/180)
            ede1 <- mad(i.coor[tp,2],center=cor[2],na.rm=T)
            ede2 <- mad(i.coor[bt,2],center=cor[2],na.rm=T)
            n1 <- length(tp)
            n2 <- length(bt)
            snr2 <- median(data$data$SNR[bt])
            snr1 <- median(data$data$SNR[tp])
            flu2 <- median(data$data$FLUX[bt])
            flu1 <- median(data$data$FLUX[tp])

          } else {
            ln2 <- ln
            cor <- apply(i.coor[ii,],2,median,na.rm=T)
            sd2 <- apply(i.coor[ii,],2,mad,na.rm=T)*c(cos(cor[2]*pi/180),1)
            tp <- intersect(ii,i1)
            bt <- intersect(ii,i2)           
            era1 <- mad(i.coor[tp,1],center=cor[1],na.rm=T)*cos(cor[2]*pi/180)
            era2 <- mad(i.coor[bt,1],center=cor[1],na.rm=T)*cos(cor[2]*pi/180)
            ede1 <- mad(i.coor[tp,2],center=cor[2],na.rm=T)
            ede2 <- mad(i.coor[bt,2],center=cor[2],na.rm=T)
            n1 <- length(tp)
            n2 <- length(bt)
            snr2 <- median(data$data$SNR[bt])
            snr1 <- median(data$data$SNR[tp])
            flu2 <- median(data$data$FLUX[bt])
            flu1 <- median(data$data$FLUX[tp])
  
            if ((era2*2 < sd2[1] | ede2*2 < sd2[2]) & n2 > 4) {
              cor <- apply(i.coor[bt,],2,median,na.rm=T)
              sd2 <- apply(i.coor[bt,],2,mad,na.rm=T)*c(cos(cor[2]*pi/180),1)
              ln2 <- length(bt)
            }
            if ((era1*2 < sd2[1] | ede1*2 < sd2[2]) & n1 > 4) {
              cor <- apply(i.coor[tp,],2,median,na.rm=T)
              sd2 <- apply(i.coor[tp,],2,mad,na.rm=T)*c(cos(cor[2]*pi/180),1)
              ln2 <- length(tp)
            }

          }
        c(cor,sd2,ln2,era1,ede1,n1,flu1,snr1,era2,ede2,n2,flu2,snr2)
        }
      
     }
  registerDoMC(1)
  mens[nde,] <- mens1  
 
return(mens)
}

mucalc <- function(data,year,weight=FALSE) {

mu <- array(dim=c(length(data[,1,1]),4),dimnames=list(star_id=attr(data,"dimnames")[[1]],data=c('mu.ra','mu.dec','mu.ra.sig','mu.dec.sig')))

ln <- length(data[1,,1])

if (weight) {
  for (i in 1:length(data[,1,1])) {
    if (length(na.omit(data[i,,3])) > 1 & length(na.omit(data[i,,4])) > 1 & !is.na(data[i,ln,3]) & !is.na(data[i,ln,4])) {
      #rr <-   lm(data[i,,1] ~ year,weights=(1/data[i,,3]*data[i,,5])^2) 
      rr <-   lm(data[i,,1] ~ year,weights=(1*data[i,,5])) 
      mu[i,1] <- (rr$coefficients[2])*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180)
      mu[i,3] <- (coef(summary(rr))[2, 2])*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180)
      if(is.na(mu[i,3])) {mu[i,3] <- sqrt(sum(data[i,,3]^2,na.rm=T))*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180) / diff(year[!is.na(data[i,,1])])}
      #tt <- lm(data[i,,2] ~ year,weights=(1/data[i,,4]*data[i,,5])^2)
      tt <- lm(data[i,,2] ~ year,weights=(1*data[i,,5]))
      mu[i,2] <- (tt$coefficients[2])*60^2 * 1000
      mu[i,4] <- coef(summary(tt))[2, 2]*60^2 * 1000
      if(is.na(mu[i,4])) {mu[i,4] <- sqrt(sum(data[i,,4]^2,na.rm=T))*60^2 * 1000 / diff(year[!is.na(data[i,,2])])}
    }
  }
} else {
  for (i in 1:length(data[,1,1])) {
    if (length(na.omit(data[i,,3])) > 1 & length(na.omit(data[i,,4])) > 1 & !is.na(data[i,ln,3]) & !is.na(data[i,ln,4])) {
      rr <-   lm(data[i,,1] ~ year)    
      mu[i,1] <- (rr$coefficients[2])*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180)
      mu[i,3] <- (coef(summary(rr))[2, 2])*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180)
      if(is.na(mu[i,3])) {mu[i,3] <- sqrt(sum(data[i,,3]^2,na.rm=T))*60^2 * 1000 * cos(mean(data[i,,2],na.rm=T)*pi/180) / diff(year[!is.na(data[i,,1])])}
      tt <- lm(data[i,,2] ~ year)
      mu[i,2] <- (tt$coefficients[2])*60^2 * 1000
      mu[i,4] <- coef(summary(tt))[2, 2]*60^2 * 1000
      if(is.na(mu[i,4])) {mu[i,4] <- sqrt(sum(data[i,,4]^2,na.rm=T))*60^2 * 1000 / diff(year[!is.na(data[i,,2])])}
    }
  }
}
return(mu)
}


correct_set <- function(data1,data2,set,sub) {
  dra <- data2[sub,1]-data1[sub,1]
  ddec <- data2[sub,2]-data1[sub,2]

data1[set,1] <- data1[set,1] + median(dra,na.rm=T)
data1[set,2] <- data1[set,2] + median(ddec,na.rm=T)
return(data1)
}

cut_error <- function(data1,data2,set,err=2) {
dra <- data2[set,1]-data1[set,1]
ddec <- data2[set,2]-data1[set,2]
bd <- which(abs(dra - mean(dra,na.rm=T)) > 2*sd(dra,na.rm=T) | abs(ddec - mean(ddec,na.rm=T)) > 2*sd(ddec,na.rm=T) )
t1 <- setdiff(set,set[bd])
return(t1)
}

select_bright <- function(data1,data2,ra,dec,members,ERROR=75,SNR=30) {
if (!is.na(ERROR)) {
  t1 <- which(!is.na(data1[,1]) & !is.na(data2[,1]) & !members & (data1$ch1_SNR > SNR | data1$ch2_SNR > SNR) & (data2$ch1_SNR > SNR | data2$ch2_SNR > SNR) & data1$e_DEC*3600*1000 < ERROR & data2$e_DEC*3600*1000 < ERROR & data1$e_RA*3600*1000*cos(data2$DEC*pi/180) < ERROR & data2$e_RA*3600*1000*cos(data2$DEC*pi/180) < ERROR & data1$DEC > dec[1] & data1$DEC < dec[2] & data1$RA > ra[1] & data1$RA < ra[2])
} else {
 t1 <- which(data1$DEC > dec[1] & data1$DEC < dec[2] & data1$RA > ra[1] & data1$RA < ra[2])
}

return(t1)
}

align <- function(data1,data2,ra,dec,exclude,ERROR=75,SNR=30) {
  t1 <- select_bright(data1,data2,ra,dec,exclude,ERROR,SNR)
  t2 <- select_bright(data1,data2,ra,dec,NA,NA,NA)
  t1 <- cut_error(data1,data2,t1)
  data1 <- correct_set(data1,data2,t2,t1)
return(data1)
}

CD.solver2 <- function(data,realz=NA,cores=4,SNR_cut=20,star_id=NA,cda=NA) {
data(wa_pars1,wa_pars2,ca_pix1)
  if (is.na(star_id[1])) {
    star_id <- data$star_id
  }

  if (length(cda) == 2) {
    sdata <- cda[[2]]
    cda <- cda[[1]]
  } else {
    sdata <- data
  }

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)

  if (mean(data$image_key$HMJD,na.rm=T) > 54976) {WARM <- TRUE} else {WARM <- FALSE}

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,sdata,1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,sdata,2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,sdata,1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,sdata,2)   
    }
  }

#number of stars:
	nst <- length(data$index)
#number of frames
	nfr <- length(data$image_key$image_id)
#median SNRs
  nde2 <- which(data$n_detections != 0)
  snr <- rep(NA,length(data$index))
  for (j in 1:length(nde2)) {
    snr[nde2[j]] <- median(data$data$SNR[data$index[[nde2[j]]]],na.rm=T)
  }

	
  nde2 <- which(snr > SNR_cut)
  nde2 <- intersect(nde2,match(star_id,data$star_id))
  star_id <- data$star_id

  realz2 <- matrix(NA,nrow=length(data$index),ncol=2)
  realz2[nde2,] <- as.matrix(realz[nde2,1:2])
  realz <- realz2
	
    if (is.na(cda[1])) {
        cda <- matrix(0, nrow = nfr, ncol = 3)
        cda[, 1:2] <- as.matrix(data[[2]][, 2:3])
        cda[, 3] <- atan2(-data[[2]][,5],data[[2]][,7])-pi
    }
	
#intial guess at coordinates:
	if (length(realz) == 1) {
		i.coor1 <- calc_all(pix1,pars1,pix2,pars2,cda,data,nde2,FALSE) 
	} else {
		i.coor1 <- realz
	}

#align to the realz.....
	#res <- single_image_resid(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id)  
	cda2 <- imageCD_solver(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id,TRUE) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}
	#res2 <- single_image_resid(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id) 


  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  nro <- length(data[[1]]$image_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda2[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],sdata,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda2[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],sdata,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  mens <- matrix(NA,nrow=length(data$index),ncol=2)
  for (j in 1:length(nde2)) {
    if (length(data$index[[nde2[j]]]) > 1) {
      mens[nde2[j],] <- apply(rresd[data$index[[nde2[j]]],],2,mean,na.rm=T)
    } else {
      mens[nde2[j],] <- rresd[data$index[[nde2[j]]],]
    }
  }

  s1 <- RMS(mens[,1])*2.5
  s2 <- RMS(mens[,1])*2.5

	
	k1 <- which(abs(mens[,1]) > s1 | abs(mens[,2]) > s2)

  realz[k1,] <- NA	
  nde2 <- setdiff(nde2,k1)
  
	cda3 <- imageCD_solver(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda3) == 4) { cda3 <- t(cda3)}
  nums <- cda3[,4]
  cda3 <- cda3[,1:3]
  if (length(cda3) == 3) { cda3 <- t(cda3)}
	#res3 <- single_image_resid(data,realz,cda3,pix1,pars1,pix2,pars2,cores,star_id) 

    i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda3[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda3[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	cda4 <- imageCD_solver(data,realz,cda3,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda4) == 4) { cda4 <- t(cda4)}
  nums <- cda4[,4]
  cda4 <- cda4[,1:3]
  if (length(cda4) == 3) { cda4 <- t(cda4)}
	#res4 <- single_image_resid(data,realz,cda4,pix1,pars1,pix2,pars2,cores,star_id) 

 i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda4[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda4[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

  cda5 <- imageCD_solver(data,realz,cda4,pix1,pars1,pix2,pars2,cores,star_id)
  if (length(cda5) == 4) { cda5 <- t(cda5)}
  nums <- cda5[,4]
  cda5 <- cda5[,1:3]
  if (length(cda5) == 3) { cda5 <- t(cda5)}
	res5 <- single_image_resid(data,realz,cda5,pix1,pars1,pix2,pars2,cores,star_id) 
		
	return(cbind(cda5,nums,res5))
}



CD.solver3 <- function(data,realz=NA,cores=4,SNR_cut=20,star_id=NA,cda=NA) {
data(wa_pars1,wa_pars2,ca_pix1)
  if (is.na(star_id[1])) {
    star_id <- data$star_id
  }

datau <- data

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  c1 <- which(data$image_key$image_id > 0)
  c2 <- which(data$image_key$image_id < 0)

  if (mean(data$image_key$HMJD,na.rm=T) > 54976) {WARM <- TRUE} else {WARM <- FALSE}

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(1,1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(1,1),2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(1,1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(1,1),2)   
    }
  }

#number of stars:
	nst <- length(data$index)
#number of frames
	nfr <- length(data$image_key$image_id)
#median SNRs
  nde2 <- which(data$n_detections != 0)
  snr <- rep(NA,length(data$index))
  for (j in 1:length(nde2)) {
    snr[nde2[j]] <- median(data$data$SNR[data$index[[nde2[j]]]],na.rm=T)
  }

	
  nde2 <- which(snr > SNR_cut)
  nde2 <- intersect(nde2,match(star_id,data$star_id))
  star_id <- data$star_id

  realz2 <- matrix(NA,nrow=length(data$index),ncol=2)
  realz2[nde2,] <- as.matrix(realz[nde2,1:2])
  realz <- realz2
	
    if (is.na(cda[1])) {
        cda <- matrix(0, nrow = nfr, ncol = 3)
        cda[, 1:2] <- as.matrix(data[[2]][, 2:3])
        cda[, 3] <- atan2(-data[[2]][,5],data[[2]][,7])-pi
    }
	
#intial guess at coordinates:
	if (length(realz) == 1) {
		i.coor1 <- calc_all(pix1,pars1,pix2,pars2,cda,data,nde2,FALSE) 
	} else {
		i.coor1 <- realz
	}

#align to the realz..... without any scale change....
	#res <- single_image_resid(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id)  
	cda2 <- imageCD_solver(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id,TRUE) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}
	#res2 <- single_image_resid(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id) 


data1 <- list()
data1$data <- data$data[i1,]
data1$image_key <- data$image_key[c1,]
data1$star_id <- data$star_id
data2 <- list()
data2$data <- data$data[i2,]
data2$image_key <- data$image_key[c2,]
data2$star_id <- data$star_id

cdd1 <- cda2[c1,]
cdd2 <- cda2[c2,]

corex1 <- realz[match(data1$data$star_id,data$star_id),1:2]
corex2 <- realz[match(data2$data$star_id,data$star_id),1:2]

xy1 <- ad2xyang(corex1[,1],corex1[,2],cdd1,data1)
ress1 <- xy1-data1$data[,3:4]

xy2 <- ad2xyang(corex2[,1],corex2[,2],cdd2,data2)
ress2 <- xy2-data2$data[,3:4]

kk <- which(!is.na(ress1[,1]) & !is.na(data1$data$x))
xfac <- lm(xy1[kk,1]~0+data1$data$x[kk])$coefficients
yfac <- lm(xy1[kk,2]~0+data1$data$y[kk])$coefficients

kk <- which(!is.na(ress2[,1]) & !is.na(data2$data$x))
xfac2 <- lm(xy2[kk,1]~0+data2$data$x[kk])$coefficients
yfac2 <- lm(xy2[kk,2]~0+data2$data$y[kk])$coefficients

f1 <- mean(c(xfac,yfac))
f2 <- mean(c(xfac2,yfac2))

data <- datau

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(f1,f1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(f2,f2),2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(f1,f1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(f2,f2),2)   
    }
  }

	cda2 <- imageCD_solver(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id,TRUE) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  nro <- length(data[[1]]$image_id)

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda2[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda2[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  mens <- matrix(NA,nrow=length(data$index),ncol=2)
  for (j in 1:length(nde2)) {
    if (length(data$index[[nde2[j]]]) > 1) {
      mens[nde2[j],] <- apply(rresd[data$index[[nde2[j]]],],2,mean,na.rm=T)
    } else {
      mens[nde2[j],] <- rresd[data$index[[nde2[j]]],]
    }
  }

  s1 <- RMS(mens[,1])*2
  s2 <- RMS(mens[,1])*2

	
	k1 <- which(abs(mens[,1]) > s1 | abs(mens[,2]) > s2)

  realz[k1,] <- NA	
  nde2 <- setdiff(nde2,k1)
  
	cda3 <- imageCD_solver(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda3) == 4) { cda3 <- t(cda3)}
  nums <- cda3[,4]
  cda3 <- cda3[,1:3]
  if (length(cda3) == 3) { cda3 <- t(cda3)}
	#res3 <- single_image_resid(data,realz,cda3,pix1,pars1,pix2,pars2,cores,star_id) 

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda3[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda3[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	cda4 <- imageCD_solver(data,realz,cda3,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda4) == 4) { cda4 <- t(cda4)}
  nums <- cda4[,4]
  cda4 <- cda4[,1:3]
  if (length(cda4) == 3) { cda4 <- t(cda4)}
	#res4 <- single_image_resid(data,realz,cda4,pix1,pars1,pix2,pars2,cores,star_id) 

 i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda4[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda4[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

  cda5 <- imageCD_solver(data,realz,cda4,pix1,pars1,pix2,pars2,cores,star_id)
  if (length(cda5) == 4) { cda5 <- t(cda5)}
  nums <- cda5[,4]
  cda5 <- cda5[,1:3]
  if (length(cda5) == 3) { cda5 <- t(cda5)}
	res5 <- single_image_resid(data,realz,cda5,pix1,pars1,pix2,pars2,cores,star_id) 


	return(list(cda=cbind(cda5,nums,res5),scale=c(f1,f2)))
}

ad2xyang <- function(ra,dec,cds,gen_set) {
##
 
  ra1 <- ra*pi/180; dec1 <- dec*pi/180
  cdec1 <- cds[match(gen_set[[1]]$image_id,gen_set[[2]]$image_id),2]*pi/180
  cra1 <- cds[match(gen_set[[1]]$image_id,gen_set[[2]]$image_id),1]*pi/180
  angs <- cds[match(gen_set[[1]]$image_id,gen_set[[2]]$image_id),3]

  phi <- pi + atan2(-cos(dec1)*sin(ra1-cra1),sin(dec1)*cos(cdec1)-cos(dec1)*sin(cdec1)*cos(ra1-cra1))
  theta <- asin(sin(dec1)*sin(cdec1)+cos(dec1)*cos(cdec1)*cos(ra1-cra1))
  rtheta <- 1/tan(theta)
  xsi <- rtheta*sin(phi)*180/pi
  eta <- -rtheta*cos(phi)*180/pi
 
  x1 <- cos(angs)*xsi + sin(angs)*eta
  y1 <- -sin(angs)*xsi + cos(angs)*eta

  return(cbind(x1,y1))
}

CD.solver4 <- function(data,realz=NA,cores=4,SNR_cut=20,star_id=NA,cda=NA) {
data(wa_pars1,wa_pars2,ca_pix1)
  if (is.na(star_id[1])) {
    star_id <- data$star_id
  }

datau <- data

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  c1 <- which(data$image_key$image_id > 0)
  c2 <- which(data$image_key$image_id < 0)

  if (mean(data$image_key$HMJD,na.rm=T) > 54976) {WARM <- TRUE} else {WARM <- FALSE}

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(1,1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(1,1),2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(1,1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(1,1),2)   
    }
  }

#number of stars:
	nst <- length(data$index)
#number of frames
	nfr <- length(data$image_key$image_id)
#median SNRs
  nde2 <- which(data$n_detections != 0)
  snr <- rep(NA,length(data$index))
  for (j in 1:length(nde2)) {
    snr[nde2[j]] <- median(data$data$SNR[data$index[[nde2[j]]]],na.rm=T)
  }

	
  nde2 <- which(snr > SNR_cut)
  nde2 <- intersect(nde2,match(star_id,data$star_id))
  star_id <- data$star_id

    if (is.na(cda[1])) {
        cda <- matrix(0, nrow = nfr, ncol = 3)
        cda[, 1:2] <- as.matrix(data[[2]][, 2:3])
        cda[, 3] <- atan2(-data[[2]][,5],data[[2]][,7])-pi
    }
	
#intial guess at coordinates:
	if (length(realz) == 1) {
		i.coor1 <- calc_all(pix1,pars1,pix2,pars2,cda,data,nde2,FALSE) 
	} else {
    realz2 <- matrix(NA,nrow=length(data$index),ncol=2)
    realz2[nde2,] <- as.matrix(realz[nde2,1:2])
    realz <- realz2
		i.coor1 <- realz
	}

#align to the realz..... without any scale change....
	#res <- single_image_resid(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id)  
	cda2 <- imageCD_solver(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id,TRUE) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}
	#res2 <- single_image_resid(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id) 


data1 <- list()
data1$data <- data$data[i1,]
data1$image_key <- data$image_key[c1,]
data1$star_id <- data$star_id
data2 <- list()
data2$data <- data$data[i2,]
data2$image_key <- data$image_key[c2,]
data2$star_id <- data$star_id

cdd1 <- cda2[c1,]
cdd2 <- cda2[c2,]

corex1 <- realz[match(data1$data$star_id,data$star_id),1:2]
corex2 <- realz[match(data2$data$star_id,data$star_id),1:2]

xy1 <- ad2xyang(corex1[,1],corex1[,2],cdd1,data1)
ress1 <- xy1-data1$data[,3:4]

xy2 <- ad2xyang(corex2[,1],corex2[,2],cdd2,data2)
ress2 <- xy2-data2$data[,3:4]

kk <- which(!is.na(ress1[,1]) & !is.na(data1$data$x))
xfac <- lm(xy1[kk,1]~0+data1$data$x[kk])$coefficients
yfac <- lm(xy1[kk,2]~0+data1$data$y[kk])$coefficients

kk <- which(!is.na(ress2[,1]) & !is.na(data2$data$x))
xfac2 <- lm(xy2[kk,1]~0+data2$data$x[kk])$coefficients
yfac2 <- lm(xy2[kk,2]~0+data2$data$y[kk])$coefficients

f1 <- mean(c(xfac,yfac))
f2 <- mean(c(xfac2,yfac2))

data <- datau

  if (WARM) {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(f1,f1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(f2,f2),2) 
    }
  } else {
    pix1 <- NA
    pix2 <- NA
    pars1 <- NA
    pars2 <- NA
    if (length(i1) > 0) {
      data$data[i1,3:4] <- pix_correc(data$data$x[i1],data$data$y[i1],ca_pix1)
      data$data[i1,3:4] <- dis_correc(data$data$x[i1],data$data$y[i1],wa_pars1,c(f1,f1),1)
    }
    if (length(i2) > 0) {
      data$data[i2,3:4] <- dis_correc(data$data$x[i2],data$data$y[i2],wa_pars2,c(f2,f2),2)   
    }
  }

	cda2 <- imageCD_solver(data,realz,cda,pix1,pars1,pix2,pars2,cores,star_id,TRUE) 
  if (length(cda2) == 4) { cda2 <- t(cda2)}
  nums <- cda2[,4]
  cda2 <- cda2[,1:3]
  if (length(cda2) == 3) { cda2 <- t(cda2)}
  #res2 <- single_image_resid(data,realz,cda2,pix1,pars1,pix2,pars2,cores,star_id)

  i1 <- which(data[[1]]$image_id > 0)
  i2 <- which(data[[1]]$image_id < 0)
  nro <- length(data[[1]]$image_id)

#new
  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda2[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],c(f1,f1),1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda2[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],c(f2,f2),2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2


  rre <- i.coor1[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])
  s2 <- RMS(rresd[,2])

	k1 <- which(abs(rresd[,1]) > 6*s1 | abs(rresd[,2]) > 6*s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor2 <- calc_all(pix1,pars1,pix2,pars2,cda2,data,data$star_id[nde2],TRUE,cores)

	cda3 <- imageCD_solver(data,i.coor2,cda2,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda3) == 4) { cda3 <- t(cda3)}
  nums <- cda3[,4]
  cda3 <- cda3[,1:3]
  if (length(cda3) == 3) { cda3 <- t(cda3)}
  	#res3 <- single_image_resid(data,i.coor2,cda3,pix1,pars1,pix2,pars2,cores,star_id) 

  i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda3[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda3[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor3 <- calc_all(pix1,pars1,pix2,pars2,cda2,data,data$star_id[nde2],TRUE,cores)
	cda4 <- imageCD_solver(data,i.coor3,cda3,pix1,pars1,pix2,pars2,cores,star_id) 
  if (length(cda4) == 4) { cda4 <- t(cda4)}
  nums <- cda4[,4]
  cda4 <- cda4[,1:3]
  if (length(cda4) == 3) { cda4 <- t(cda4)}
	#res4 <- single_image_resid(data,realz,cda4,pix1,pars1,pix2,pars2,cores,star_id) 

 i.coor.1 <- coor.calc(pix=pix1,pars=pars1,crsa=cda4[match(data[[1]]$image_id[i1],data[[2]]$image_id),],u=data[[1]]$x[i1],v=data[[1]]$y[i1],data,1)
  i.coor.2 <- coor.calc(pix=pix2,pars=pars2,crsa=cda4[match(data[[1]]$image_id[i2],data[[2]]$image_id),],u=data[[1]]$x[i2],v=data[[1]]$y[i2],data,2)

  i.coor <- matrix(NA,nrow=nro,ncol=2)
  i.coor[i1,] <- i.coor.1; i.coor[i2,] <- i.coor.2

  rre <- realz[match(data[[1]]$star_id,star_id),1:2]
  rresd <- ((rre[,1:2] - i.coor)*cbind(cos(rre[,2]*pi/180),rep(1,length(i.coor[,1])))*3600)

  s1 <- RMS(rresd[,1])*2
  s2 <- RMS(rresd[,2])*2

	k1 <- which(abs(rresd[,1]) > s1 | abs(rresd[,2]) > s2)
  
  data$data$x[k1] <- NA
  data$data$y[k1] <- NA

	i.coor4 <- calc_all(pix1,pars1,pix2,pars2,cda3,data,data$star_id[nde2],TRUE,cores)
  cda5 <- imageCD_solver(data,i.coor4,cda4,pix1,pars1,pix2,pars2,cores,star_id)
  if (length(cda5) == 4) { cda5 <- t(cda5)}
  nums <- cda5[,4]
  cda5 <- cda5[,1:3]
  if (length(cda5) == 3) { cda5 <- t(cda5)}
	res5 <- single_image_resid(data,i.coor4,cda5,pix1,pars1,pix2,pars2,cores,star_id) 


	return(list(cda=cbind(cda5,nums,res5),scale=c(f1,f2)))
}

