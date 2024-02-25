function a = a_irT(i,r,T)
% Return Helmholtz function (J/kg) of species i given density (kg/m3) and temp. (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
a = R_i(i)*T*( a0_idt(i,d,t) + ar_idt(i,d,t) );