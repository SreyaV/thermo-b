function u = u_irT(i,r,T)
% Return internal energy (J/kg) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
u = R_i(i)*T*( t*(a0t_idt(i,d,t) + art_idt(i,d,t)) );