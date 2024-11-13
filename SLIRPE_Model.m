function [dydt] = SLIRPE_Model(t,y,params, T, epsilon)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% params = [farm.beta, farm.muL, farm.muI, farm.k, T(1), Tmin, Tmax];


% index params
beta = params(1);
muL = params(2);
muI = params(3);
k = params(4);
T = params(5);
Tmin = params(6);
Tmax = params(7);
%fprintf('K=%d',k)
Tunit=T-Tmin;                  %change later
   % if Tmin < T && T<Tmax
   %      Tunit=T-Tmin;
   %      fprintf('Tunit=%d',Tunit)
   % elseif T<Tmin || T>Tmax
   %          Tunit=0;
   %          %fprintf('T_unit= %d', Tunit)
   % end

S = y(1);
L = y(2);
I = y(3);
R = y(4);
P = y(5);


%%
dE = epsilon;
dP=  k*Tunit+dE;
dS=  -((beta*S*I)+dE) +dP;
dL= ((beta*S*I)+dE) -(L*(muL));
dI=  (L*muL) -(I*muI);
dR=  (I*muI);
dydt = [dS; dL; dI; dR; dP;dE];
%%  %stuff for finding dF
% Kf is a constnt
%TunitStar find from integral below eq 5 pdf
% M=sqrt(U^2+V^2)
f_sporeReleaseNewInfect=@(M) 1/(1+(Mmax/nu)*exp(1-M))
F_growth= k_f * TunitStar* I

dF = F_growth- F_growth*f_sporeReleaseNewInfect;

end
