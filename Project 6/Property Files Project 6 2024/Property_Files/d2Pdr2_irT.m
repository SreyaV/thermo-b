function d2Pdr2 = d2Pdr2_irT(i,r,T)
% Return the second derivative of P wrt r (Pa/kg/m3)^2 for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);

% Make an ideal solution.

f = 2*d*ard_idt(i,d,t) + ...
    4*d*d*ardd_idt(i,d,t) + ...
    d*d*d*arddd_idt(i,d,t);

% Scale.

d2Pdr2 = R_i(i)*T*f/r;