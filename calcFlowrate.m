function [ q, Re ] = calcFlowrate( L, rho, P_wh, mu, D, P_wf, e )
%CALCFLOWRATE calculates single-phase flow-rate in a pipe
%   considering turbulent flow and in field units the flowrate is obtained
%   via Chen's correlation (1979). Returns flow-rate and the Reynolds
%   number.

gc = 32.17;

syms u f
Re = 1488.8*rho*u*D/mu;

eqnf = 1/sqrt(f) == -4*log10(e/3.7065 - 5.0452/Re*log10((e^1.1098)/2.8257...
    + (7.149/Re)^.8981)); 
f = solve(eqnf, f);

eqn = (P_wf - P_wh)*144 == rho*L + 2*f*rho*u^2*L/gc/D;
u = solve(eqn, u);

q = double(u)*pi*D^2/4/0.000065; % bbl/day
Re = 1488.8*rho*double(u)*D/mu;

end