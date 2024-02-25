function Cv = Cv_crT(c,r,T)
% Return the heat capacity (J/kg-K) for r (kg/m3) and T (K);
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.
f = 0;
for i=1:1:length(c)
    f = f + c(i)*(...
        -(Tcrit_i(i)/T)^2*a0tt_idt(i,r/rcrit_i(i),Tcrit_i(i)/T)...
        - t^2*artt_idt(i,d,t)...
        );
end

% Add the excess.
f = f - t^2*aett_cdt(c,d,t);

% Scale.
Cv = R_c(c)*f;
