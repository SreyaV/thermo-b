function s = s_crT(c,r,T)
% Return entropy (J/kg-K) for given c, r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 0;
for i=1:1:length(c)
    if (c(i) ~= 0)
        f = f + c(i)*((Tcrit_i(i)/T)*a0t_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) + t*art_idt(i,d,t)...
            - (a0_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) + ar_idt(i,d,t)) - log(c(i)));
    end
end


% Add the excess.

f = f + t*aet_cdt(c,d,t) - ae_cdt(c,d,t);

% Scale.

s = R_c(c)*f;