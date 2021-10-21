
library(gmm)
library(AER)
library(stargazer)

us = read.csv('./Output/3_US_1.csv')
cpi_dum = ifelse(us$cpi[2:101] < us$cpi[1:100], 1, 0)
proxy_dum = ifelse(us$proxy < 0, 1, 0)[2:101]
gdp.pprox = (us$rgdpprox[2:101] - us$rgdpprox[1:100]) / us$rgdpprox[1:100] * 100
gdp.cpi = (us$rgdpcpi[2:101] - us$rgdpcpi[1:100]) / us$rgdpcpi[1:100] * 100

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

MODELS = function(yt,xt,zt) {
  # Baseline(OLS, CPI)
  m1 = lm(yt~xt)
  
  # Independent(OLS, CPI, Proxy)
  x0z0 = (1-xt)*(1-zt)
  x1z0 = xt*(1-zt)
  x0z1 = (1-xt)*zt
  x1z1 = xt*zt
  m2 = lm(yt~x0z1+x1z0+x1z1)
  
  # Independent(GMM, CPI, Proxy)
  t0 = c(m2$coefficients[1],m2$coefficients[4],p=0.5,ex=0.15,ey=0.15,vx=0.15,vy=0.15)
  x = as.matrix(cbind(yt,x0z0,x1z0,x0z1,x1z1))
  m3 = gmm(g=g1, x=x, t0=t0, type="iterative")
  
  # Independent(IV, CPI, Proxy)
  m4 = ivreg(yt~xt|zt)
  
  # Independent(IV, CPI, Proxy)
  t0 = c(m2$coefficients[1],m2$coefficients[4],p=0.5,e11=0.15,e10=0.15,v10=0.15,v00=0.15)
  m5 = gmm(g2, x, t0, type="iterative")
  return(list(m1,m2,m3,m4,m5))
}

m.cpi = MODELS(gdp.cpi,cpi_dum,proxy_dum)
m.pprox = MODELS(gdp.pprox,cpi_dum,proxy_dum)

for(i in 1:5) {
  stargazer(m.cpi[[i]], m.pprox[[i]], type = 'text')
}

# Simulation
set.seed(4274)
n = 1000
cpi = 2.8*rnorm(n)
cpi_dum = ifelse(cpi<0, 1, 0)

beta=-1
alpha=3
gdp = alpha + beta*cpi_dum + 2*rnorm(n)
ngdp = gdp + cpi

errx = -0.3 + 3*rnorm(n)
errz = -0.3 + 3*rnorm(n)
errg = -0.3 + 3*rnorm(n)

cpix = cpi + errx
cpiz = cpi + errz
cpig = cpi + errg

gdpx = ngdp - cpix
gdpz = ngdp - cpiz
gdpg = gdp + errg

cpi_x_dum = ifelse(cpix<0, 1, 0)
cpi_z_dum = ifelse(cpiz<0, 1, 0)
cpi_g_dum = ifelse(cpig<0, 1, 0)

m.well = MODELS(gdp, cpi_x_dum, cpi_z_dum)
m.err1 = MODELS(gdpg, cpi_x_dum, cpi_z_dum)
m.err2 = MODELS(gdpx, cpi_x_dum, cpi_z_dum)

for(i in 1:5) {
  stargazer(m.well[[i]], m.err1[[i]], m.err2[[i]], type = 'text')
}



