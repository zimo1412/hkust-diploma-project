# United State
us = read.csv("./Data/RawDataUS.csv", skip = 1)
us = data.frame(us)
us = us[us$Year >= 1799 & us$Year < 1900, ]

d.19 = read.csv('./Output/1_Historical.csv')
d.19 = data.frame(d.19)
d.19 = d.19[d.19$date >= 1799 & d.19$date < 1900,]


#Create error-ridden real GDP for the US
NGDP = us$RGDP*us$GDPD
RGDPPROXY = NGDP/d.19$realProxy
RGDPCPI = NGDP/us$CPI

library(mFilter)
FiltGDP = hpfilter(log(us$RGDP[!is.na(us$RGDP)]), freq = 100, type=c("lambda"))
FiltIP = hpfilter(log(us$IPROD[!is.na(us$IPROD)]),freq = 100, type=c("lambda"))
GAPIP = FiltIP$cycle*100
GAPGDP = FiltGDP$cycle*100

# Export historical data for use in Stata
us.data = cbind(us$Year, d.19$Proxy, us$RGDP, RGDPPROXY, RGDPCPI, 
                us$GDPD, us$RGDPPC, GAPIP, GAPGDP, us$CPI, us$CPIF, us$BANK, 
                us$IPROD, us$IPFood, us$IPTextiles, us$IPLeather, us$IPChemicals, 
                us$IPMachinery, us$IPWood, us$IPMetals, us$SHARE, us$M2, us$M3)
colnames(us.data) = c("year", "proxy", "rgdp", "rgdpprox", "rgdpcpi", 
                      "gdpd", "rgdppc", "gapip", "gap", "cpi", "cpif", "bank", 
                      "iprod", "ipfood",  "iptex", "ipleat", "ipchem", 
                      "ipmach", "ipwood", "ipmet", "share", "m2", "m3")
write.csv(us.data, "./Output/3_US_1.csv", row.names=F)


# United Kingdom
uk = read.csv("./Data/RawDataUK.csv", skip = 1)
uk = data.frame(uk)

# Construct proxy
Warable = uk$Warable[1]
Wpasture = uk$Wpasture[1]
Wwood = uk$Wwood[1]
uk = uk[uk$Year >= 1830 & uk$Year < 1900, ]

Parable = uk$Parable
Ppasture = uk$Ppasture
Pwood = uk$Pwood
# Calculate linked series
library(zoo)
Parable = na.approx(Parable, na.rm = F)
Ppasture = na.approx(Ppasture, na.rm = F)
Pwood = na.approx(Pwood, na.rm = F)

# Compute proxy
P = cbind(Parable, Ppasture, Pwood)
W = rbind(Warable, Wpasture, Wwood)
w.m = function(p) { weighted.mean(p,W,na.rm=T) }
PROX = apply(P, 1, w.m)

# Export historical data for use in Stata
uk.data = cbind(uk$Year, uk$CPI, uk$WPI, PROX, uk$GDEFL, uk$CDEFL, uk$GDP, 
               uk$GDPPC, uk$PROD, uk$CONS, uk$INV, uk$GCONS, uk$EXP, uk$IMP, 
               uk$URATE, uk$STX, uk$HOUSE)
colnames(uk.data) = c("year", "cpi", "wpi", "proxy", "gdefl", "cdefl", "gdp", 
                      "gdppc", "prod", "cons", "inv", "gcons", "exp", "imp", 
                      "urate", "stx", "house")
write.csv(uk.data, "./Output/3_UK_1.csv", row.names=F)
