function d3Pdr2dT = d3Pdr2dT_irT(i,r,T)
% Return the third derivative of P wrt r and T (Pa/(kg/m3)^2/K) for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);

% Make an ideal solution.

f = (2*d*ard_idt(i,d,t) + ...
     4*d*d*ardd_idt(i,d,t) + ...
     d*d*d*arddd_idt(i,d,t))...
  - (Tcrit_i(i)/T)*(2*d*ardt_idt(i,d,t) + ...
     4*d*d*arddt_idt(i,d,t) + ...
     d*d*d*ardddt_idt(i,d,t));

% Scale.

d3Pdr2dT = R_i(i)*f/r;