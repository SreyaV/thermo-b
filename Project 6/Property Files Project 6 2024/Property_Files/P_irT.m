function P = P_irT(i,r,T)
% Return the pressure (Pa) given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i R_i

t  = Tcrit_i(i)/T;
d  = r/rcrit_i(i);

P  = r*R_i(i)*T*( 1 + d*ard_idt(i,d,t) );