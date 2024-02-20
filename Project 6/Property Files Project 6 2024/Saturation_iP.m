function [Tsat rf rg] = Saturation_iP(ispecies,P)
% Return the saturation Temperature (K) and density (kg/m3) for any species i 
% at any given P (Pa).

Setup_Props_i;

% First make a guess using Dühring plot

% Critical and triple points
Pcrit = Pcrit_i(ispecies);
Tcrit = Tcrit_i(ispecies);
Ptrip = Ptrip_i(ispecies);
Ttrip = Ttrip_i(ispecies);

% ln(P) vs -1/T is affine : ln(P) = a * (-1/T) + b
a = (log(Pcrit) - log(Ptrip)) / (-1/Tcrit -  (-1/Ttrip));
b = log(Pcrit) - a * ( -1/Tcrit );

% Guess is (-1/T) = (log(P) - b) / a
inv_T  = (log(P)-b) / a;
Tguess = -1/inv_T;


% Then use saturation_iT with bisection method since
% Pressure against Temperature is increasing

T_upper = min(Tguess + 2,0.99*Tcrit);
T_lower = max(Tguess - 2,Ttrip);
tolerance = 1e-3;

while T_upper - T_lower > tolerance

    T = (T_upper + T_lower)/2;
    [Psat rfsat rgsat] = Saturation_iT(ispecies,T);
   
    if Psat < P
        T_lower = T;
    else
        T_upper = T;
    end
end

Tsat = (T_upper + T_lower)/2;
rf   = rfsat;
rg   = rgsat;
