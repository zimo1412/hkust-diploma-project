us = read.csv('./Output/3_US_1.csv')
cpi_dum = ifelse(us$cpi[2:101] < us$cpi[1:100], 1, 0)
proxy_dum = ifelse(us$proxy < 0, 1, 0)[2:101]
d.iprod = (us$iprod[2:101] - us$iprod[1:100]) / us$iprod[1:100] * 100

# Baseline(OLS, CPI)
m1 = lm(d.iprod~cpi_dum)
summary(m1)

# Independent(OLS, CPI, Proxy)
c0p0 = (1-cpi_dum)*(1-proxy_dum)
c0p1 = (1-cpi_dum)*proxy_dum
c1p0 = cpi_dum*(1-proxy_dum)
c1p1 = cpi_dum*proxy_dum

m2 = lm(d.iprod~c0p1+c1p0+c1p1)
summary(m2)

# Independent(GMM, CPI, Proxy)
library(gmm)
g1 = function(theta,data) {
  a = theta[1]; b = theta[2]
  p = theta[3]
  ex = theta[4]; ez = theta[5]
  vx = theta[6]; vz = theta[7]
  yt = data[,1]
  x0z0 = data[,2]
  x1z0 = data[,3]
  x0z1 = data[,4]
  x1z1 = data[,5]
  
  y.x1z1 = x1z1 * (yt - (a + b*p*(1-ex)*(1-ez)/(p*(1-ex)*(1-ez)+vx*vz*(1-p))))
  y.x0z1 = x0z1 * (yt - (a + b*p*ex*(1-ez)/(p*ex*(1-ez)+(1-vx)*vz*(1-p))))
  y.x1z0 = x1z0 * (yt - (a + b*p*(1-ex)*ez/(p*(1-ex)*ez+vx*(1-vz)*(1-p))))
  y.x0z0 = x0z0 * (yt - (a + b*p*ex*ez/(p*ex*ez+(1-vx)*(1-vz)*(1-p))))
  
  Px1z1 = x1z1 - (1-ex)*(1-ez)*p + vx*vz*(1-p)
  Px1z0 = x1z0 - ex*(1-ez)*p + (1-vx)*vz*(1-p)
  Px0z1 = x0z1 - (1-ex)*ez*p + vx*(1-vz)*(1-p)
  
  f = cbind(y.x1z1,y.x1z0,y.x0z1,y.x0z0,Px1z1,Px1z0,Px0z1)
  return(f)
}

data = as.matrix(cbind(d.iprod,c0p0,c1p0,c0p1,c1p1))
set.seed(1)
t0 = c(a=m2$coefficients[1], b=m2$coefficients[4], p=0.5, 
       ex=0.15, ez=0.15, vx=0.15, vz=0.15)

m3 = gmm(g=g1, x=data, t0=t0, type="twoStep")
summary(m3)

# Independent(IV, CPI, Proxy)
library(AER)
m4 = ivreg(d.iprod~cpi_dum|proxy_dum)
summary(m4)

# Dependent(GMM, CPI, Proxy)
g2 = function(theta,data) {
  a = theta[1]; b = theta[2]
  p = theta[3]
  e11 = theta[4]; e10 = theta[5]; 
  fix.e00 = 0.15; fix.v11 = 0.15
  v10 = theta[6]; v00 = theta[7]
  yt = data[,1]
  x0z0 = data[,2]
  x1z0 = data[,3]
  x0z1 = data[,4]
  x1z1 = data[,5]
  
  y.x1z1 = x1z1 * (yt - (a + b*p*e11/(p*e11+fix.v11*(1-p))))
  y.x0z1 = x0z1 * (yt - (a + b*p*(1-e11-e10-fix.e00)/(p*(1-e11-e10-fix.e00)+
                                                        (1-fix.v11-v10-v00)*(1-p))))
  y.x1z0 = x1z0 * (yt - (a + b*p*e10/(p*e10+v10*(1-p))))
  y.x0z0 = x0z0 * (yt - (a + b*p*fix.e00/(p*fix.e00+v00*(1-p))))
  
  Px1z1 = x1z1 - e11*p + fix.v11*(1-p)
  Px1z0 = x1z0 - (1-e11-e10-fix.e00)*p + (1-fix.v11-v10-v00)*(1-p)
  Px0z1 = x0z1 - e10*p + v10*(1-p)
  
  f = cbind(y.x1z1,y.x1z0,y.x0z1,y.x0z0,Px1z1,Px1z0,Px0z1)
  return(f)
}

data = as.matrix(cbind(d.iprod,c0p0,c1p0,c0p1,c1p1))
set.seed(1)
t0 = c(a=m2$coefficients[1], b=m2$coefficients[4], p=0.5, 
       e11=0.15, e10=0.15, v10=0.15, 
       v00=0.15)

m5 = gmm(g=g2, x=data, t0=t0, type="iterative")
summary(m5)
