inj
ob

d[1]-ob*pve[1,4]
d[2]-ob*pve[2,4]
d[3]-ob*pve[3,4]


vaxall=function(par,d, ob, pve){
  alvax=par[1:3]
  alvax <- alvax/sum(alvax)
  nd=d-ob*alvax*pve[,4]
  return(sum(nd))
}
alvax=c(0.33, 0.33, 0.33)
 pars=c(alvax,sum(alvax))
out <- optim(par=pars, # initial guess
   fn=vaxall,d=d,ob=ob,pve=pve,
   method="L-BFGS-B",lower=c(0,0,0,.9),upper=c(.9,.9,.9,1))

P_init=alvax
my_ui=rbind(c(0,1), c(0, 1), c(0,1))
z = constrOptim(P_init,vaxall,NULL,ui=my_ui, ci=my_ci)

######## solve as transportation problem


# Import lpSolve package
library(lpSolve)

# Set efficacy matrix 
eff <- t(pve)*100

# Set unequality/equality signs for suppliers
row.signs <- rep("<=", 4)

inj=5000 # # of injectable vaccines acquired
obh=16000
# Set right hand side coefficients for suppliers
row.rhs <- c(sum(d), sum(d), sum(d), sum(d))

# Set unequality/equality signs for customers
col.signs <- rep("<=", 3)

# Set right hand side coefficients for # dogs
col.rhs <- c(d[1], d[2], d[3])*0.7

# Final value (z)
lp.transport(eff, "min", row.signs, row.rhs, col.signs, col.rhs)

# Variables final values
lp.transport(eff, "max", row.signs, row.rhs, col.signs, col.rhs)$solution

#### try to add more constrainst using regular lp function
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
injc=c(1,1,1,0,1,1,1,0,1,1,1,0)
obhc=c(0,0,0,1,0,0,0,1,0,0,0,1)
constr=rbind(constr, injc, obhc)
rhs <- c(d[1], d[2], d[3], inj, inj, inj, ob, inj, ob)
constr.dir <- c(rep("<=",m ), rep("<=", n+2))
#rhs <- c(d[1], d[2], d[3], inj, inj, inj, ob)
#constr.dir <- c(rep("<=",m ), rep("<=", n))

prod.trans <- lp ("max", objective.in=obj, 
                  constr, constr.dir, rhs)
soln=matrix(prod.trans$solution, nrow=4)

colSums(soln)/d
