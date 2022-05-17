---
title: "optimize rabies campaign: take 1"
author: "Kim Cuddington"
date: "07/05/2022"
output: 
  html_document: 
    keep_md: yes
---



# Optimizing rabies vaccine programs for dogs

Our goal is to use standard optimization tools to find the most cost-effective rabies vaccination strategy. This doc is take 1: get an optimal vaccination plan for a given number of vaccines and dogs.

We're using information from:

[Wallace, R. M., Undurraga, E. A., Gibson, A., Boone, J., Pieracci, E. G., 
Gamble, L., & Blanton, J. D. (2019). Estimating the effectiveness of 
vaccine programs in dog populations. Epidemiology & Infection, 147.](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E0D25E5DFB352731121EDBA5DEFEF7EE/S0950268819001158a.pdf/estimating-the-effectiveness-of-vaccine-programs-in-dog-populations.pdf)

online version [https://rabiestaskforce.com/toolkit/vaxplan/](https://rabiestaskforce.com/toolkit/vaxplan/)

# Dog population structure

The dog population is thought of as consisting of always confined dogs (C), sometimes confined dogs (SC), and never confined (NC). The status is independent of ownership, and the population of free-roaming dogs is therefore SC+NC.

If we know the total dog population and the percent in each category, we can also estimate numbers this way.

To start, we're just using the default values from the online tool.


```r
# Dt=total dog population 
# pC[i] = proportion of dogs in each category
# Categories are C always confined, SC sometimes confined, NC never confined
# and sum to one

Dt=30000
C=0.28
NC=0.48
SC=1-(C+NC)
pC=c(C,SC,NC)
d=Dt*pC
```


# Accessibility

Wallace et al (2019) have a category called accessibility,
which seems to relate to how easy it is to administer the vaccine, using the delivery method. These methods are central point (CP), door-to-door (DD), capture, vaccinate and release (CVR), and oral bait handout (ORV). The first three methods use an injectable vaccine.



```r
# vaccine efficacy estimates
# pve = probability of vaccinating a dog in category C using method i
# methods are CP - central point, DD-door to door, 
# CVR-capture, vaccinate, release, ORV-oral bait
# values from Wallance (assume ORV is more effective with never confined)

CP=c(.95,.80,.05)
DD=c(.95,.80,.05)
CVR=c(.05,.95,.80)
ORV=c(.05,.95,.95)

pve=data.frame(CP,DD,CVR,ORV)
vaxcats=colnames(pve)
dogcats=c("C", "SC", "NC")
rownames(pve)=dogcats
```

# Vaccine Solutions

At this point Wallace t al. (2019) describe an "iterative set of equations to 
define the number of vaccinated dogs in each category by vaccine method, using 
number of of vaccines. A pre-defined rank order is defined 
to determine which category would be most accessible for a given method
and allocated vaccines applied to that category first, before applying remainder 
to the next category, until all vaccines used. Excess is waste"

For example, the spreadsheet gives the following ranking


```r
#vaccine ranking
C_V=c(1,1,3,3)
SC_V=c(2,2,1,2)
NC_V=c(3,3,2,1)

pvr=t(data.frame(C_V,SC_V,NC_V))
rownames(pve)=dogcats
colnames(pve)=vaxcats
pvr
```

```
     [,1] [,2] [,3] [,4]
C_V     1    1    3    3
SC_V    2    2    1    2
NC_V    3    3    2    1
```

I'm not sure why this is "predefined", since we can certainly turn the
given efficacy into rank. Although I also note a small discrepancy
in the ranking of SC given by the spreadsheet, so maybe there is another factor here.


```r
# turn efficacy into rank
rank.pve=pve
ad=cut(as.matrix(rank(-pve)), 3, labels=c(1:3))
rank.pve[] <- as.numeric(ad)
rank.pve
```

```
   CP DD CVR ORV
C   1  1   3   3
SC  2  2   1   1
NC  3  3   2   1
```

The online tool then starts with the idea that a given number of vaccines 
have been obtained, and they need to be allocated using an iterative process
embodied in their spreadsheet according to a predefined vaccine strategy.


```r
#vaccines obtained
#St = total number of vaccine doses obtained
#pm = proportion of doses allocated to each vaccination method 
inj=20400
ob=3600
```

Rather than examining a pre-determined vaccination strategy that may or may not be successful, it seems much simpler to simply allocate the vaccines to different methodss given the probability of success for the various dog categories using linear programming. 

However, I do note that when this is done, there is no particular priority given to
any dog category, so that we do need to define the end goal (i.e., 70% coverage in each category, or only in some high priority category). Otherwise, some optimal solutions end up with only 1 or 2 categories of dogs being vaccinated when there are not enough vaccines to go around. 

Other constraints include:
1. not exceeding the total number of vaccines, and 2. not exceeding the total
number of dogs in each category. 

Please note for reference that this is a "transportation problem" with additional constraints.


```r
#### optimize with constraints for a given number of 
#### vaccines, using regular lp function
library(lpSolve)
obj=c(t(as.matrix(pve)))*100
m <- 3
n <- 4
constr <- matrix (0 , n +m , n*m )
for(i in 1:m){ 
  for(j in 1:n){ 
    constr[i, n*(i-1) + j] <- 1
    constr[m+j, n*(i-1) + j] <- 1
  }
}
# this array will ensure we cannot exceed the number of dogs in 
# a category, and that a given the number of vaccines in a category
# cannot exceed the total

injc=c(1,1,1,0,1,1,1,0,1,1,1,0) # all categories summed cannot exceed inj vaccines
obhc=c(0,0,0,1,0,0,0,1,0,0,0,1) #all categories summed cannot exceed obh vaccines
cover1=c(c(1,1,1,1,0,0,0,0,0,0,0,0)) #required coverage C
cover2=c(c(0,0,0,0,1,1,1,1,0,0,0,0)) #required coverage SC
cover3=c(c(0,0,0,0,0,0,0,0,1,1,1,1)) #required coverage NC

constr=rbind(constr, injc, obhc, cover1, cover2, cover3)
rhs <- c(d[1], d[2], d[3], inj, inj, inj, 
         ob, inj, ob, 0.7*d[1],0.7*d[2],0.7*d[3])
constr.dir <- c(rep("<=",m ), rep("<=", n+2),
                rep(">=",3))


prod.trans <- lp ("max", objective.in=obj, 
                  constr, constr.dir, rhs)
soln=(matrix(prod.trans$solution, nrow=4))
solntab=(rbind(soln,colSums(soln)/d ))
colnames(solntab)=dogcats
rownames(solntab)=c(vaxcats, "%vax")
knitr::kable(solntab, digits=2, caption="Table: Optimal vaccine allocation strategy")
```



Table: Table: Optimal vaccine allocation strategy

|     |    C|      SC|     NC|
|:----|----:|-------:|------:|
|CP   | 8400|    0.00|    0.0|
|DD   |    0|    0.00|    0.0|
|CVR  |    0| 5520.00| 6480.0|
|ORV  |    0|    0.00| 3600.0|
|%vax |    1|    0.77|    0.7|

This solution gives us good coverage, but I bet the expense of CVR makes it untenable. Next up... incorporate uncertainty.


## Uncertainty


Confidence intervals for these predictions are not easy to generate analytically (*ref). Instead, we will use a Monte Carlo method to quantify uncertainty. To do this we will use the beta distribution. This is a continuous probability distribution that models random variables with values falling inside a finite interval, and we can use it to model the success rates of the vaccination methods, since these have both an upper and lower bound for possible values. 

The beta distribution has two shape parameters, α and β. Both parameters must be positive values. The distribution is particularly flexible at modeling different curves within the interval, including symmetrical, left and right-skewed, U and inverted U shapes, and straight lines. Using its tight relationship to the binomial distribution, with 10 trials and 7 successes we would have: 
    α = 7 + 1 = 8,
    β = 10 – 7 + 1 = 4

For example, 


```r
## Visualization, including limit cases:
pl.beta <- function(a,b, asp = if(isLim) 1, ylim = if(isLim) c(0,1.1)) {
  if(isLim <- a == 0 || b == 0 || a == Inf || b == Inf) {
    eps <- 1e-10
    x <- c(0, eps, (1:7)/16, 1/2+c(-eps,0,eps), (9:15)/16, 1-eps, 1)
  } else {
    x <- seq(0, 1, length.out = 1025)
  }
  fx <- cbind(dbeta(x, a,b), pbeta(x, a,b), qbeta(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
  matplot(x, f, ylab="", type="l", ylim=ylim, asp=asp,
          main = sprintf("[dpq]beta(x, a=%g, b=%g)", a,b))
  abline(0,1,     col="gray", lty=3)
  abline(h = 0:1, col="gray", lty=3)
  legend("top", paste0(c("d","p","q"), "beta(x, a,b)"),
         col=1:3, lty=1:3, bty = "n")
  invisible(cbind(x, fx))
}

aval=7+1
bval=10-7+1

aval
```

```
[1] 8
```

```r
bval
```

```
[1] 4
```

```r
pl.beta(aval, bval)
```

![](notesoptimrabies_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


Then we use this distribution with many simulations


```r
# function to generate random vaccine efficacy matrix
rbpve<-function(x){
aval=x*10+1
bval=10-x*10+1
y=rbeta(1, aval, bval)
return(y)
}

# function to calculate optimal allocation with random matrix draw
# many times
randvax<-function(pve,d) {
covmat=matrix(NA, nrow=100, ncol=3)
mclist=list()
plist=list()
for (i in 1:100){
  p=apply(pve, c(1,2), rbpve)
  plist[[i]]=p
  obj=c(t(as.matrix(p)))*100
  prod.trans <- lp ("max", objective.in=obj,constr, constr.dir, rhs)

  soln=(matrix(prod.trans$solution, nrow=4))
mclist[[i]]=soln
covmat[i,]=colSums(soln)/d
}
return(list(mclist, covmat, plist))
}
multilist=randvax(pve,d)

#function to tabulate the resulting simulations

tabprep<-function(flist) {
mclist=multilist[1]
plist=multilist[3]
#  Make a 3D array from list of matrices
arr <- array( unlist(mclist) , c(4,3,100) )
parr<-array( unlist(plist) , c(3,4,100) )

#  Get summaries of third dimension
mdarr=apply( arr , 1:2 , median)
mxarr=apply( arr , 1:2 , max)
mnarr=apply( arr , 1:2 , min)

mdp=round(apply( parr , 1:2 , median), 2)
mxp=round(apply( parr , 1:2 , max), 2)
mnp=round(apply( parr , 1:2 , min), 2)
pmat=matrix(paste(mdp, "(", mnp, "-", mxp, ")"), nrow=3,ncol=4)
row.names(pmat)=dogcats
omat=t(matrix(paste(mdarr, "(", mnarr, "-", mxarr, ")"), nrow=3,ncol=4, byrow=T))
covmat=matrix(unlist(multilist[2]),  ncol=3)
cmax=apply(covmat, 2, max)
cmin=apply(covmat, 2, min)
cmed=apply(covmat, 2, median)
crow=paste(cmed, "(", cmin, "-", cmax, ")")
omat=rbind(omat, crow)
colnames(omat)=dogcats
rownames(omat)=c(vaxcats, "%vax")


return(list(pmat, omat))

}
tablist=tabprep(multilist)
pmat=tablist[1]

#Print the results
knitr::kable(pmat, digits=2, caption="Table: Median and range of vaccine efficacy values", col.names=vaxcats)
```



<table class="kable_wrapper">
<caption>Table: Median and range of vaccine efficacy values</caption>
<tbody>
  <tr>
   <td> 

|   |CP                   |DD                   |CVR                  |ORV                  |
|:--|:--------------------|:--------------------|:--------------------|:--------------------|
|C  |0.9 ( 0.62 - 0.99 )  |0.9 ( 0.54 - 1 )     |0.1 ( 0 - 0.65 )     |0.11 ( 0.01 - 0.57 ) |
|SC |0.76 ( 0.32 - 0.98 ) |0.79 ( 0.47 - 0.96 ) |0.91 ( 0.61 - 1 )    |0.9 ( 0.56 - 1 )     |
|NC |0.12 ( 0.01 - 0.4 )  |0.1 ( 0 - 0.36 )     |0.77 ( 0.38 - 0.94 ) |0.9 ( 0.57 - 1 )     |

 </td>
  </tr>
</tbody>
</table>

```r
omat=tablist[2]
knitr::kable(omat, digits=2, caption="Table: Optimal vaccine allocation strategy (median, min - max)", col.names=dogcats)
```



<table class="kable_wrapper">
<caption>Table: Optimal vaccine allocation strategy (median, min - max)</caption>
<tbody>
  <tr>
   <td> 

|     |C                 |SC                            |NC                              |
|:----|:-----------------|:-----------------------------|:-------------------------------|
|CP   |6720 ( 0 - 8400 ) |0 ( 0 - 7200 )                |0 ( 0 - 0 )                     |
|DD   |0 ( 0 - 8400 )    |0 ( 0 - 7200 )                |0 ( 0 - 0 )                     |
|CVR  |0 ( 0 - 0 )       |5520 ( 0 - 7200 )             |6480 ( 6480 - 13080 )           |
|ORV  |0 ( 0 - 0 )       |0 ( 0 - 3600 )                |3600 ( 0 - 3600 )               |
|%vax |1 ( 0.7 - 1 )     |0.766666666666667 ( 0.7 - 1 ) |0.7 ( 0.7 - 0.908333333333333 ) |

 </td>
  </tr>
</tbody>
</table>

It's not clear to me that the best way to portray the range of outcomes....



We may also have uncertainty in the number of dogs, and/or the number of dogs in each category.

We can use a beta distribution for proportions and a standard normal distribution for the total numbers.


```r
#Let's assume a 20% CV for total dog numbers, and no randomness in proportions to start

sddogs=.2*Dt

for (j in 1:10){
  nDt=rnorm(1, Dt, sddogs)
print(nDt)
d=nDt*pC
multilist=randvax(pve,d)
}
```

```
[1] 22717.71
[1] 31967.26
[1] 29866.3
[1] 38124.6
[1] 28351.84
[1] 31287.45
[1] 24769.49
[1] 29828.51
[1] 32558.02
[1] 34658.81
```

```r
# not satisfactory since we have to run 100*# reps of dogs
# could just change function manually, or...


#idea to change function def on the fly to incorporate dog # variability
# first save the definition as a list of string
##newDef <- deparse(vis.gam)

# then identify the line to be changed using regular expressions
# (see ?regexp)
##iLine <- grep("gray\\(seq\\(",initDef)

# replace the line by what you want
##newDef[iLine] <- "            pal <- gray(seq(0.9, 0.1, length = nCol))"

# and define a new function by parsing and evaluating the 
# new definition
##vis.gam2 <- eval(parse(text=newDef))

## Next up
#change the constraints from number of vax units to $$
```
