function h = h0_crT(c,r,T)
% Return the ideal gas enthalpy (J/kg) for given c, r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

% Make an ideal solution.
f = 1;
for i=1:1:length(c)
    if(c(i) ~= 0)
        f = f + c(i)*(Tcrit_i(i)/T)*a0t_idt(i,r/rcrit_i(i),Tcrit_i(i)/T);
    end
end

% Scale.
h = R_c(c)*T*f;
