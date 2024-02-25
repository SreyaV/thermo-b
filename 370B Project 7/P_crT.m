function P = P_crT(c,r,T)
% Return the pressure (Pa) given r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 1;
for i=1:1:length(c)
    f = f + c(i)*d*ard_idt(i,d,t);
end

% Add the excess.

f = f + d*aed_cdt(c,d,t);

% Scale.

P = r*R_c(c)*T*f;