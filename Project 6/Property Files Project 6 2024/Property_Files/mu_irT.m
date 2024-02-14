function mu = mu_irT(i,r,T)
% Return the chemical potential (J/kmol) for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i M_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
g = R_i(i)*T*( 1 + a0_idt(i,d,t) + ar_idt(i,d,t)+ d*ard_idt(i,d,t) );
mu = g*M_i(i);