function s0 = s0_irT(i,r,T)
% Return the ideal gas entropy (J/kg-K) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
s0 = R_i(i)*( t*a0t_idt(i,d,t) - a0_idt(i,d,t) );