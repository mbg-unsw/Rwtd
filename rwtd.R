# R translation of the Stata wtdttt package
# https://github.com/mbg-unsw/Rwtd
#
# (c) 2023 Malcolm Gillies <malcolm.gillies@unsw.edu.au>
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>. 

# This code is a direct translation of portions from the Stata module
# wtdttt, available at:
# https://econpapers.repec.org/software/bocbocode/s458265.htm
# by Katrine Bødkergaard Nielsen (kani@ph.au.dk) and
#    Henrik Stovring (stovring@ph.au.dk) 

library(data.table)
library(bbmle)

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
	-sum(log(plogis(logitp) * pnorm(-(log(obstime) - mu)/exp(lnsigma)) /
	    exp(mu + exp(2*lnsigma)/2) + plogis(-logitp) / delta))
    }

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

## wtdttt rx1time, disttype(lnorm) start(1jan2014) end(31dec2014)

# library(haven)
#
# x <- read_dta("ref/wtddat_dates.dta")
#
# x.w <- wtdttt(x$rx1time, as.Date('2014-01-01'), as.Date('2014-12-31'),
#	reverse=F)
# 
# summary(x.w)
# 
# x <- read_dta("ref/ranwtddat_discdates.dta")
# 
# x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'),
#	as.Date('2014-12-31'), reverse=F)
# 
# summary(x.w)
# 
# x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'), as.Date('2014-12-31'),
#     reverse=T)
# 
# summary(x.w)
