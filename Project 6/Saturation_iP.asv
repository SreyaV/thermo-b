function [Tsat rf rg] = Saturation_iP(ispecies,P)
% Return the saturation Temperature (K) and density (kg/m3) for any species i 
% at any given P (Pa).


% First make a guess using Dühring plot

% Critical and triple points
Pcrit = Pcrit_i(ispecies);
Tcrit = Tcrit_i(ispecies);
Ptrip = Ptrip_i(ispecies);
Ttrip = Ttrip_i(ispecies);

% ln(P) vs -1/T is affine : ln(P) = a * (-1/T) + b
a = (ln(Pcrit) - ln(Ptrip)) / (-1/Tcrit -  (-1/Ttrip));
b = ln(Pcrit) - a * ( -1/Tcrit );

% Guess is (-1/T) = (ln(P) - b) / a
inv_T  = (ln(P)-b) / a;
Tguess = -1/inv_T;


% Then use saturation_iT with bisection method since
% Pressure against Temperature is increasing

T_upper = Tguess + 1;
T_lower = Tguess - 1;
tolerance = 1e-3;

while T_upper - T_lower > tolerance

    T = (T_upper + T_lower)/2;
    [Psat rfsat rgsat] = Saturation_iT(ispecies,T);

    if Psat < P
        T_upper = T;
    else
        T_lower = T;
    end
end

Tsat = (T_upper + T_lower)/2;
rf   = rfsat;
rg   = rgsat;
