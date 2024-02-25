function g = g_crT(c,r,T)
% Return the Gibbs function (J/kg) for r (kg/m3), c and T (K).
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 1;
for i=1:1:length(c)
    f = f + c(i)*(a0_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) + ar_idt(i,d,t) + d*ard_idt(i,d,t));
end

% Add the excess.

f = f + ae_cdt(c,d,t)+ d*aed_cdt(c,d,t);

% Scale.

g = R_c(c)*T*f;