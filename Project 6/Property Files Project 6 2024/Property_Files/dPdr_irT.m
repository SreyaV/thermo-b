function dPdr = dPdr_irT(i,r,T)
% Return the derivative of P wrt r (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
dPdr = R_i(i)*T*( 1 + 2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t) );