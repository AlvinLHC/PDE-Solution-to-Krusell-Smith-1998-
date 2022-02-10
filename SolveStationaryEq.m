function residual = SolveStationaryEq(r) 

global alpha delta lmean 
steady = stationary_eq(r);
CapitalSupply = steady.A;
CapitalDemand = (alpha/(r+delta))^(1/(1-alpha))*lmean;
residual = CapitalSupply - CapitalDemand;
