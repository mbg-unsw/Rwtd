# Rwtd
**n.b. this code has been superseded by the [Rwtdttt](https://github.com/mbg-unsw/Rwtdttt) package**

Rough R translation of functions from the Stata wtdttt package

(c) 2023 Malcolm Gillies <malcolm.gillies@unsw.edu.au>

This code is a direct translation of portions from the Stata module
wtdttt, available at:
https://econpapers.repec.org/software/bocbocode/s458265.htm
by Katrine Bødkergaard Nielsen (kani@ph.au.dk) and
   Henrik Støvring (stovring@ph.au.dk) 

These functions estimate parameters of the ordinary and reverse Waiting Time
Distribution (WTD) by maximum likelihood, with application to modelling the
recurrence of drug dispensing. Please see the Stata package documentation
for references.

As with wtdttt, this code is made available under the terms of the GPL v3 (or later).

*Currently, only the lognormal distribution is implemented.*

Example data in Stata format, taken directly from the wtdttt package,
is included for demonstration purposes

```R
library(haven)

x <- read_dta("ref/wtddat_dates.dta")

x.w <- wtdttt(x$rx1time, as.Date('2014-01-01'), as.Date('2014-12-31'),
      reverse=F)

summary(x.w)

x <- read_dta("ref/ranwtddat_discdates.dta")

x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'),
      as.Date('2014-12-31'), reverse=F)

summary(x.w)

x.w <- ranwtdttt(x$rxdate, x$pid, as.Date('2014-01-01'), as.Date('2014-12-31'),
    reverse=T)

summary(x.w)
```
