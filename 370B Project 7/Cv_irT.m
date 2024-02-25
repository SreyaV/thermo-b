function Cv = Cv_irT(i,r,T)
% Return the heat capacity (J/kg-K) for r (kg/m3) and T (K);
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
Cv = R_i(i)*( -t*t*(a0tt_idt(i,d,t) + artt_idt(i,d,t)) );