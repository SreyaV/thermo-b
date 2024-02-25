function s = s0_crT(c,r,T)
% Return the ideal gas entropy (J/kg-K) for given c, r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

global Tcrit_i rcrit_i

% Make an ideal solution.
f = 0;
for i=1:1:length(c)
    if(c(i) ~= 0)
        f = f + c(i)*...
            (...
            (Tcrit_i(i)/T)*a0t_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) ...
            - a0_idt(i,r/rcrit_i(i),Tcrit_i(i)/T) - log(c(i))...
            );
    end
end

% Scale.
s = R_c(c)*f;
