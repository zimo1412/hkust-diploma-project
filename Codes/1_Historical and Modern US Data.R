cpi = read.csv('./Data/CompositeCPIMeasuringWorth.csv')
proxy = read.csv('./Data/RawDataProxy.csv')


library(lubridate)
library(xts)
library(tsbox)
library(seasonal)
load(file="./Data/BLSData.RData")

calcIndex <- function(series, weights, baseY) {
  # Useful function to calculate weighted mean of indexed series
  
  series <- ts_index(series, baseY)*100
  
  Index <- sapply(seq_len(nrow(series)), function(i){
    weighted.mean(as.matrix(series[i,]), weights, na.rm = TRUE)
  })
  
  Index <- xts(Index, order.by = as.Date(index(series)))
  
  return(Index)
}

# F1: Flour
F1w = cbind(9.7, 0.04)
startNew = as.Date(index(f12[1]))
startOld = as.Date(index(f11[1]))
F1 = ts_chain(ts_span(f10, NULL, startOld), f11)
F1 = ts_chain(ts_span(F1, NULL, startNew), f12)

# F2: Rice, pasta, cornmeal
F2w = cbind(1.3, 0.119)
startNew = as.Date(index(f22[1]))
startOld = as.Date(index(f21[1]))
f21b = ts_chain(ts_span(f20, startOld), f21)
F2 = ts_chain(ts_span(f21b, NULL, startNew), f22)

# F3: Beef and veal
F3w = cbind(7.9, 0.436)

# F4: Pork
F4w = cbind(2.8, 0.3)
F4 = ts_chain(f41, f42)

# F5: Mutton
F5w = cbind(1.5, 0.227)
startNew = as.Date(index(f51[1]))
F5 = ts_chain(ts_span(f50, NULL, startNew), f51)

# F6: Fish and seafood
F6w = cbind(1.4, 0.25)

# F7: Fresh whole milk
F7w = cbind(3, 0.207)

# F8: Butter
F8w = cbind(5.9, 0.0305)

# F9: Cheese
F9w = cbind(0.5, 0.242)
startNew = as.Date(index(f91[1]))
F9 = ts_chain(ts_span(f90, NULL, startNew), f91)

# F10: Fresh vegetables
F10w = cbind(7.0, 0.478)

# F11: Fresh fruit
# Seasonally adjust potatoes to estimate relative variance of growth rates
prt = final(seas(ts_ts(prt)))
pws = final(seas(ts_ts(pws)))
pwsSd = sd(ts_pc(pws), na.rm=TRUE)
prtSd = sd(ts_pc(prt), na.rm=TRUE)
pwsMn = mean(ts_pc(pws), na.rm=TRUE)
prtMn = mean(ts_pc(prt), na.rm=TRUE)
pwsAdj = (ts_pc(pws)-pwsMn)/pwsSd*prtSd+pwsMn/pwsSd*prtSd
pwsAdjLev = ts_compound(pwsAdj)
fwsSd = sd(ts_pc(fws), na.rm =TRUE)
frtSd = sd(ts_pc(frt), na.rm =TRUE)
fwsMn = mean(ts_pc(fws), na.rm =TRUE)
fwsAdj = (ts_pc(fws)-fwsMn)/pwsSd*prtSd+fwsMn/pwsSd*prtSd
fwsAdjLev = ts_compound(fwsAdj)

F11w = cbind(2.5, 0.555)
F11 = fwsAdjLev


# F12: Eggs
F12w = cbind(2.3, 0.097)

# F13: Other beverages (Tea)
F13w = cbind(1.3, 0.093)
F13 = ts_chain(f131, f132)

# F14: Coffee
F14w = cbind(4.0, 0.168)
F14 = ts_chain(f141, f142)

# F15: Margarine
F15w = cbind(1.5, 0.0305)

# F16: Sugar and sweets
F16w = cbind(4.8, 0.282)

# C1: Women's apparel
C1w = cbind(3.4, 1.044)
C1 = ts_chain(ts_span(c11, NULL, ts_summary(c12)$start), c12)

# C2: Men's apparel
C2w = cbind(3.6, 0.581)
C2 = ts_chain(ts_span(c21, NULL, ts_summary(c22)$start), c22)

# C3: Footwear
C3w = cbind(4.0, 0.671)
C3 = ts_chain(ts_span(c31, NULL, ts_summary(c32)$start), c32)

# H1: Personal care products (Soap)
H1w = cbind(0.3, 0.704)
H1 = ts_bind(h11, h14)
H1 = ts_xts(na.approx(ts_ts(H1)))

# H2: Housekeeping supplies (Starch)
H2w = cbind(0.2, 0.835)
startNew = as.Date(index(h22[1]))
startOld = as.Date(index(h21[1]))
valOld = h21[startOld]
h21b = ts_chain(ts_span(h20, NULL, startOld), h21)
valNew = h21b[startOld]
h21b = h21b/as.numeric(valNew[1])*as.numeric(valOld[1])
H2 = ts_bind(h21b, h22)
H2 = ts_xts(na.approx(ts_ts(H2)))

# H3: Textile furnishings
# ------------------------------------------------------------------------
H3w = cbind(3.5, 0.261)
H3 = ts_bind(ts_index(h31, "1997-12-01")*100, h32)
H3 = ts_xts(na.approx(ts_ts(H3)))

# E1: Energy commodities (Fuel and Light)
E1w = cbind(7, 4.679)
startNew = as.Date(index(e11[1]))
valOld = e11[startNew]
e11b = ts_chain(ts_span(e10b, NULL, startNew), e11)
valNew = e11b[startNew]
e11b = e11b/as.numeric(valNew[1])*as.numeric(valOld[1])
E1 = ts_bind(e11b, e12)
E1 = ts_xts(na.approx(ts_ts(E1)))

# O1: Newspapers
O1w = cbind(1.1, 0.067)
O1 = ts_bind(ts_index(o11, "1997-12-01")*100, o12)
O1 = ts_xts(na.approx(ts_ts(O1)))

# O2: Repair of household items (Shoe repair)
O2w = cbind(0.7, 0.104)
startNew = as.Date(index(o21[1]))
o21b = ts_chain(ts_span(o20, NULL, startNew), o21)
O2 = ts_chain(o21b, o22)

# O3: Medical care services
O3w = cbind(1.1, 6.924)
O3 = ts_chain(o31, o32)

# R1: Rent
R1w = cbind(17.7, 30.099)

# ppiF: Food
ppiFw = cbind(57.4, 13.384)

# ppiC: Clothing (Apparel)
ppiCw = cbind(11.0,3.018)

# ppiH: Household commodities (textile house furnishings)
ppiHw = cbind(4.0, 1.778)

# ppiE: Energy (fuels and related products and power)
ppiEw = cbind(7.0,	4.634)

# ppiHh: Household furnishings historical missing data (link with best match)
startNew = as.Date(index(ts_na_omit(ppiHh)[1]))
ppiHh = ts_chain(ts_span(ppiCh, NULL, startNew), ppiHh)

# CPI replication
P = ts_c(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, C1, C2, C3, H1, H2, H3, E1, O1, O2, R1)
W = rbind(F1w, F2w, F3w, F4w, F5w, F6w, F7w, F8w, F9w, F10w, F11w, F12w, F13w, F14w, F15w, F16w, C1w, C2w, C3w, H1w, H2w, H3w, E1w, O1w, O2w, R1w)
rownames(W) = colnames(P)
colnames(W) = c("Hoover", "2018")

CPIw19 = calcIndex(P, W[, 1], baseY = "2010-12-01")
CPIw21 = calcIndex(P, W[, 2], baseY = "2010-12-01")

CPIw19 = ts_span(CPIw19, "1960-01-01")
CPIw21 = ts_span(CPIw21, "1960-01-01")

# Actual CPI and PPI data from BLS
CPI = ts_span(CPI, "1960-01-01")
PPI = ts_span(PPI, "1960-01-01")

# Proxy
P  = ts_c(ppiF, ppiC, ppiH, ppiE)
Ph = ts_c(ppiFh, ppiCh, ppiHh, ppiEh)
W  = rbind(ppiFw, ppiCw, ppiHw, ppiEw)
rownames(W) = colnames(P)
colnames(W) = c("Hoover", "2018")

# Modern replications
proxyw19 = calcIndex(P, W[, 1], baseY = "2010-12-01")
proxyw21 = calcIndex(P, W[, 2], baseY = "2010-12-01")

proxyw19 = ts_span(proxyw19, "1960-01-01")
proxyw21 = ts_span(proxyw21, "1960-01-01")

# Historical proxy
proxyhw19 = calcIndex(Ph, W[, 1], baseY = "1890-01-01")

# Plots for paper in annual frequency (Modern)
CPIa = ts_pcy(ts_frequency(CPI, to="year", aggregate = "mean"))
PROXa = ts_pcy(ts_frequency(proxyw21, to="year", aggregate = "mean"))
REPa = ts_pcy(ts_frequency(CPIw21, to="year", aggregate = "mean"))
date = as.Date(index(REPa))

data.df = data.frame(date, CPIa, PROXa, REPa)
colnames(data.df) = c("date", "CPI", "Proxy", "Replication")
write.csv(data.df,'./Output/1_Modern.csv', row.names = F)

# Plots for paper in annual frequency (Historical)
CPIa = ts_pcy(cpiAllh)
PROXa = ts_pcy(proxyhw19)
date = as.Date(index(CPIa))

data.df = data.frame(date, ts_c(CPIa, PROXa, proxyhw19))
colnames(data.df) = c("date", "CPI", "Proxy","realProxy")
write.csv(data.df,'./Output/1_Historical.csv', row.names = F)

