function h0 = h0_irT(i,r,T)
% Return the ideal gas enthalpy (J/kg) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
h0 = R_i(i)*T*( 1 + t*a0t_idt(i,d,t) );