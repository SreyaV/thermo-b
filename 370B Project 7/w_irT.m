function w = w_irT(i,r,T)
% Return speed of sound (m/s) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global R_i

t  = Tcrit_i(i)/T;
d  = r/rcrit_i(i);
w2 = R_i(i)*T*( 1 + 2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t) - ...
    (1 + d*ard_idt(i,d,t) - d*t*ardt_idt(i,d,t))^2/...
    (t*t*(a0tt_idt(i,d,t) + artt_idt(i,d,t))) );
w  = sqrt(w2);