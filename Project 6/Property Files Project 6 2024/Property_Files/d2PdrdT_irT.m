function d2PdrdT = d2PdrdT_irT(i,r,T)
% Return the 2nd derivative of P wrt r and T (Pa/kg/m3/K) for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);
d2PdrdT = R_i(i)*( 1 + 2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t) )...
        - R_i(i)*(Tcrit_i(i)/T)*( 2*d*ardt_idt(i,d,t) + d*d*arddt_idt(i,d,t) );