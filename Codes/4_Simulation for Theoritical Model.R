sims = function(alpha, beta, sig, rho, c, N, s2n, dep, prefix) {
  set.seed(4274)
  plim = NULL
  for (i in 1:length(s2n)) {
    # (1) calculate variance
    sig.omega = sqrt(sig^2*(1/(1+s2n[i])))
    sig.pi = sqrt(sig^2 - sig.omega^2)
    # (2) calculate omega and pi
    omega = rnorm(N,0,sig.omega)
    pit = rnorm(N,0,sig.pi)
    pit.w = rho[1] + rho[2]*pit + omega
    # (3) calculate x and d
    dt = pit < c
    xt = pit.w < c
    # (4) calculate plim
    Pd1x0 = mean(dt[which(xt==F)]==T)
    Pd0x1 = mean(dt[which(xt==T)]==F)
    plim.a = alpha + beta*Pd1x0
    plim.b = beta*(1 - Pd0x1 - Pd1x0)
    # check if add dependent variable
    if (dep == T) {
      dpi = pit.w - pit
      Epix0 = mean(dpi[which(xt==F)])
      Epix1 = mean(dpi[which(xt==T)])
      plim.a = plim.a - Epix0
      plim.b = plim.b + Epix0 - Epix1
    }
    plim = rbind(plim, c(plim.a, plim.b))
  }
  # set colnames
  colnames(plim) = c(paste0(prefix,'_alpha'), paste0(prefix,'_beta'))
  return(plim)
}

N = 1000000
s2n = seq(0.001, 6, by=0.05)

# Test for 5 settings
plim1 = sims(alpha=1, beta=-1, sig=6, rho=c(0,1), c=0, N=N, s2n=s2n, dep=F, prefix='Baseline')
plim2 = sims(alpha=1, beta=-1, sig=6, rho=c(0,1), c=5, N=N, s2n=s2n, dep=F, prefix='Larger threshold')
plim3 = sims(alpha=1, beta=-1, sig=6, rho=c(5,1), c=0, N=N, s2n=s2n, dep=F, prefix='Larger interpret')
plim4 = sims(alpha=1, beta=-1, sig=6, rho=c(0,3), c=0, N=N, s2n=s2n, dep=F, prefix='Larger slope')
plim5 = sims(alpha=1, beta=-1, sig=6, rho=c(0,1), c=0, N=N, s2n=s2n, dep=T, prefix='Correlated')

result = cbind(s2n,plim1,plim2,plim3,plim4,plim5)
write.csv(result,'./Output/2_Simulation_1.csv', row.names = F)
