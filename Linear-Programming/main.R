rm(list = ls())
library(tidyverse)
#library(quantmod)
#library(Matrix)



#
#  Specify targets
#
T1 = c(0.30, 0.09, 0.07, 0.20, 0.08, 0.06, 0.03, 0.03, 0.02, 0.02, 0.01, 0.02, 0.07, 0.00, 0.00)
T2 = c(0.24, 0.09, 0.06, 0.11, 0.05, 0.05, 0.03, 0.03, 0.01, 0.01, 0.01, 0.03, 0.26, 0.00, 0.02)
T3 = c(0.16, 0.06, 0.04, 0.08, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.02, 0.04, 0.35, 0.04, 0.08)

#
#  Indicate which funds (by specifying indexes) 
#  that each client currently holds
#
Zachary = c(1, 4, 13, 14, 15)
Yolanda = c(3, 5, 7, 8, 9, 11, 12)

##################################################################
#   Here is where you change values to consider a different case
#
    CurrentExposures = Zachary
    Target = T1
#
##################################################################



#
#  Build constraint matrix
#


# Full investment constraint
e = as.matrix(rep(1,15))

# Primary Constraints
PChar1 = as.matrix(c(rep(1,12),rep(0,3)))
PChar3 = as.matrix(c(rep(0,7),rep(1,3),rep(0,4),1))

# Secondary Constraints
#SChar1 = as.matrix(c(rep(1,3),rep(0,3),0.55, 0.60,rep(0,7)))
#SChar2 = as.matrix(c(0,0,1,rep(0,3),0.1,1,rep(0,7)))
SChar1 = as.matrix(c(rep(1,3),rep(0,2),0.55, 0.60,rep(0,8)))
SChar2 = as.matrix(c(0,0,1,rep(0,2),0.1,1,rep(0,8)))
SChar3 = as.matrix(c(rep(0,7),rep(1,3),rep(0,5)))
SChar4 = as.matrix(c(rep(0,10),1,1,rep(0,3)))
SChar5 = as.matrix(c(rep(0,13),1,0))

SChar = t(cbind(SChar1, SChar2, SChar3, SChar4, SChar5))

#
# Full constraint matrix
#
#  Rows:
#    1. Full allocation
#    2. Primary characteristic 1
#    3. Primary characteristic 3
#    4-8. Secondary characteristic differences
#    9-23.Fund allocation differences
#    24-38. Binary allocation constraints
#    39 Sum of binary variables
#
#  Columns:
#    1-15 Allocation to each fund
#    16-20 d+ deficit of client portfolio secondary characteristics versus target
#    21-25 d- surplus of client portfolio secondary characteristics versus target
#    26-40 delta+ deficit of client portfolio fund allocations versus target
#    41-55 delta- surplus of client portfolio fund allocations versus target
#    56-70 binary allocation indicator variables
#
AMat = rbind(t(e), t(PChar1),t(PChar3), SChar, diag(15),diag(15))
AMat <- cbind(AMat,rbind(matrix(0,nrow=3,ncol = 55),
                         cbind(diag(5),-diag(5),matrix(0,nrow = 5, ncol = 45)),
                         cbind(matrix(0,nrow=15,ncol=10),diag(15),-diag(15),matrix(0,nrow=15,ncol=15)),
                         cbind(matrix(0,nrow=15,ncol=40),-diag(15))))
AMat <- rbind(AMat,cbind(matrix(0,nrow=1,ncol=55),t(e)))


## Convert AMat to sparse represenation
rownames(AMat) <- c(1:dim(AMat)[1])
colnames(AMat) <- c(1:dim(AMat)[2])
A <- as.data.frame(as.table(AMat))
colnames(A) <- c("row","col","value")
A <- A %>% filter(value != 0)


#
#  Objective coefficients
#
ObjVec <- c(rep(0,15),rep(3,10),rep(2,30),rep(0,15))

#
#  Right hand side values
#
RHS <- c(1,
         t(PChar1) %*% Target,
         t(PChar3) %*% Target, 
         SChar %*% Target,
         t(Target),
         rep(0,15),
         length(CurrentExposures)+1)

#
#  Variable upper and lower bounds
# 

lower <- c(rep(0,70))
upper <- c(rep(1,70))
lower[CurrentExposures+55] = 1

##
##
##  Let's solve the model
##
##

## solve using glpk via the API
library(glpkAPI)
## initialize model
lp<- initProbGLPK()
## model dimensions
nrows <- dim(AMat)[1]
ncols <- dim(AMat)[2]

# row upper and lower bounds
rlower <- RHS
rlower[24:38] <- -1
rlower[39] <- length(CurrentExposures)
rupper <- RHS

# maximize objective GLP_Min (maximize with GLP_MAX)
setObjDirGLPK(lp,GLP_MIN)

# tell model how many rows and columns
addRowsGLPK(lp, nrows)
addColsGLPK(lp, ncols)

# indicate variables that are integer
getColsKindGLPK(lp,1:70)
setColsKindGLPK(lp,1:ncols, c(rep(GLP_CV,55), rep(GLP_IV,15)))
getColsKindGLPK(lp,1:70)

#getColTypeGLPK(lp, 1:70)
#getColsStatGLPK(lp)
#getColsKindGLPK(lp,1:70)
#setColKindGLPK(lp,70,GLP_BV)


# add column limits
setColsBndsGLPK(lp,c(1:ncols), lower, upper)
setRowsBndsGLPK(lp,c(1:nrows),rlower,rupper)
setObjCoefsGLPK(lp,c(1:ncols),ObjVec)

# load constraint matrix coefficients
loadMatrixGLPK(lp,nrow(A), A$row, A$col,A$value)

# solve LP problem using Simplex Method
solveSimplexGLPK(lp)
solveMIPGLPK(lp)

# get retsults of MIP solution
mipStatusGLPK(lp)
status_codeGLPK(mipStatusGLPK(lp))

# report results
mipObjValGLPK(lp)
sm <- mipColsValGLPK(lp)
mipRowsValGLPK(lp)
sm[1:15]
sm[56:70]

# set integer variables to optimal values
# and resolve with simplex to get sensitivity 
# results
setColsBndsGLPK(lp,c(56:70), sm[56:70], upper)
solveSimplexGLPK(lp)
# objective function value
getObjValGLPK(lp)
# value of variables in optimal solution
getColsPrimGLPK(lp)
# status of each variable in optimal solution 1 = basic variable
getColsStatGLPK(lp)
getRowsDualGLPK(lp)
# get dual values for each row/constraint
getRowDualGLPK(lp,1)
getRowDualGLPK(lp,2)
getRowDualGLPK(lp,3)

getColsPrimGLPK(lp)[1:15]
getColsPrimGLPK(lp)[16:20]
getColsPrimGLPK(lp)[21:25]



