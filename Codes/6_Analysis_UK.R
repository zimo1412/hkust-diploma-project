
library(gmm)
library(AER)
library(stargazer)

uk = read.csv('./Output/3_UK_1.csv')
cpi_dum = ifelse(uk$cpi[2:70] < uk$cpi[1:69], 1, 0)
proxy_dum = ifelse(uk$wpi[2:70] < uk$wpi[1:69], 1, 0)
d.gdp = (uk$gdp[2:70] - uk$gdp[1:69]) / uk$gdp[1:69] * 100
d.prod = (uk$prod[2:70] - uk$prod[1:69]) / uk$prod[1:69] * 100
d.urate = uk$urate[2:70]
d.cons = (uk$cons[2:70] - uk$cons[1:69]) / uk$cons[1:69] * 100
d.inv = (uk$inv[2:70] - uk$inv[1:69]) / uk$inv[1:69] * 100

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

m.g = MODELS(d.gdp,cpi_dum,proxy_dum)
m.p = MODELS(d.prod,cpi_dum,proxy_dum)
m.u = MODELS(d.urate,cpi_dum,proxy_dum)
m.c = MODELS(d.cons,cpi_dum,proxy_dum)
m.i = MODELS(d.inv,cpi_dum,proxy_dum)

for(i in 1:5) {
  stargazer(m.g[[i]], m.p[[i]], m.u[[i]], m.c[[i]], m.i[[i]], type = 'text')
}
