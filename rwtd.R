#x <- read_dta("ref/wtddat_dates.dta")

#program define mlwtdttt_lnorm
#        version 14.0
#	args lnf logitp mu lnsigma 
#
#qui{
#        replace `lnf' = ln(invlogit(`logitp') *                   /*
#                        */ normal(-(ln($ML_y1) - `mu')/exp(`lnsigma')) /*
#                        */ / exp(`mu' + exp(2 * `lnsigma')/2)  /*
#                        */ + invlogit(-`logitp') / $wtddelta )
#      }
#end


# wtdttt rx1time, disttype(lnorm) start(1jan2014) end(31dec2014)

# first try lnorm, discrete time version only

wtdttt <- function(t, start, end, reverse) {
    # continuity correction: start - 0.5, end + 0.5
    delta <- as.double(end - start, units="days") + 1
    if(reverse)
	obstime <- 0.5 + as.double(end - t, units="days")
    else
	obstime <- 0.5 + as.double(t - start, units="days")

    # compute initial values
    ntot <- length(t)
    nonprevend <- sum(obstime > (delta * 2/3))
    prp <- 1 - 3 * nonprevend / ntot
    lpinit <- qlogis(prp)

    # inits for lognormal
    muinit <- mean(log(obstime[obstime < 0.5 * delta]))
    lnsinit <- log(sd(log(obstime[obstime < 0.5 * delta])))

    mllnorm <- function(mu, lnsigma, logitp) {
#	message('mu=', mu, '; lnsigma=', lnsigma, '; logitp=', logitp)
	-sum(log(plogis(logitp) * pnorm(-(log(obstime) - mu)/exp(lnsigma)) /
	    exp(mu + exp(2*lnsigma)/2) + plogis(-logitp) / delta))
    }

#        replace `lnf' = ln(invlogit(`logitp') *                   /*
#                        */ normal(-(ln($ML_y1) - `mu')/exp(`lnsigma')) /*
#                        */ / exp(`mu' + exp(2 * `lnsigma')/2)  /*
#                        */ + invlogit(-`logitp') / $wtddelta )
#      }

#        (logitp: `obstime' = `logitpcovar' `allcovar') ///
#        (mu: `mucovar' `allcovar') ///
#        (lnsigma: `lnsigmacovar' `allcovar') `if' `in', iterate(`niter') ///
# niter=50

    mle2(mllnorm,start=list(mu=muinit,lnsigma=lnsinit,logitp=lpinit))
}

ranwtdttt <- function(t, id, start, end, nsamp=1, reverse=T) {
    delta <- as.double(end - start, units="days")
    x <- data.table(t=t, id=id, key=c("id", "t"))
    off <- data.table(id=unique(id), key=c("id"))
    st <- as.Date(x = integer()) # empty Date vector
    for(i in 1:nsamp) {
	off[, inddate:=start + floor(runif(n=nrow(off), max=delta+1))]

	x[off, inddate := i.inddate]

	if(reverse) {
	    tmp <- x[t >= inddate - delta & t <= inddate, .SD[.N], by=id]
	    tmp[, t := t + (end - inddate)]
	} else {
	    tmp <- x[t <= inddate + delta & t >= inddate, .SD[1L], by=id]
	    tmp[, t := t - (inddate - start)]
	}
	st <- c(st, tmp$t) # XXXX FIXME
    }

    wtdttt(st, start, end, reverse)
}

# x.w <- wtdttt(x$rx1time, as.Date('2014-01-01'), as.Date('2014-12-31'),
#     reverse=F)
# 
# summary(x.w)
# 
# x <- read_dta("ref/ranwtddat_discdates.dta")
# 
# x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'), as.Date('2014-12-31'),
#     reverse=F)
# 
# summary(x.w)
# 
# x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'), as.Date('2014-12-31'),
#     reverse=T)
# 
# summary(x.w)
