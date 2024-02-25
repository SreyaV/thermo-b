function u = u_crT(c,r,T)
% Return internal energy (J/kg) for given c, r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.
f = 0;
for i=1:1:length(c)
    f = f + c(i)*((Tcrit_i(i)/T)*a0t_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) ...
        + t*art_idt(i,d,t));
end

% Add the excess.
f = f + t*aet_cdt(c,d,t);

% Scale.
u = R_c(c)*T*f;
