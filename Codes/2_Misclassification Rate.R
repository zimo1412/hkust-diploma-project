historical = read.csv('./Output/1_Historical.csv')

# Historical Statistics
CPIa = historical$CPI
CPIa.19 = CPIa[3:103]
mean(CPIa.19, na.rm=TRUE)
sd(CPIa.19, na.rm=TRUE)

PROXa = historical$Proxy
PROXa.19 = PROXa[3:103]
mean(PROXa.19, na.rm=TRUE)
sd(PROXa.19, na.rm=TRUE)

cor(CPIa.19, PROXa.19, use = "pairwise.complete.obs")

def_proxy = PROXa.19 < 0
def_cpi = CPIa.19 < 0
inf_proxy = PROXa.19 >= 0
inf_cpi = PROXa.19 >= 0
mean(def_proxy != def_cpi, na.rm = TRUE)

# Modern Misclassification Rate

MISCLASS_TABLE = function(d, x, z) {
  Px0d1 = mean(x[d==1]==0, na.rm = T)
  Px1d0 = mean(x[d==0]==1, na.rm = T)
  
  Pz0d1 = mean(z[d==1]==0, na.rm = T)
  Pz1d0 = mean(z[d==0]==1, na.rm = T)
  
  Pxz0d1 = mean(z[d==1]==0 & x[d==1]==0, na.rm = T)
  Pxz1d0 = mean(z[d==0]==1 & x[d==0]==1, na.rm = T)
  
  Table = as.matrix(cbind(Px0d1, Px1d0, Pz0d1, Pz1d0, Px0d1*Pz0d1,  Px1d0*Pz1d0, Pxz0d1, Pxz1d0))
  return(Table)
}

modern = read.csv('./Output/1_Modern.csv')
CPIa.21 = modern$CPI
PROXa.21 = modern$Proxy
REPa.21 = modern$Replication
# Zero Threshold
d = CPIa.21 < 0
x = PROXa.21 < 0
z = REPa.21 < 0
zero.thres = MISCLASS_TABLE(d,x,z)

# Higher Threshold
mCPIa = mean(CPIa.21, na.rm = TRUE)
d = CPIa.21 < mCPIa
x = PROXa.21 < mCPIa
z = REPa.21 < mCPIa
higher.thres = MISCLASS_TABLE(d,x,z)

# Save result
result = rbind(zero.thres, higher.thres)
colnames(result) = c("$x=0|d=1$","$x=1|d=0$", "$z=0|d=1$","$z=1|d=0$", "prod|d=0", "prod|d=1" , "$x=z=0|d=1$","$x=z=1|d=0$")
rownames(result) = c("1960-2017 zero threshold", "1960-2017 higher threshold")
write.csv(result, "./Output/1_Mis_Rate.csv")
