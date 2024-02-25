function s = s_irT(i,r,T)
% Return entropy (J/kg-K) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
s = R_i(i)*( t*(a0t_idt(i,d,t) + art_idt(i,d,t)) - (a0_idt(i,d,t) + ar_idt(i,d,t)) );