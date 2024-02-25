function ar = ar_crT(c,r,T)
% Return the residual Helmholtz function (J/kg) for r (kg/m3), c and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.
f = 0;
for i=1:1:length(c)
    if(c(i) ~= 0)
        f = f + c(i)*(ar_idt(i,d,t) + log(c(i)));
    end
end

% Add the excess.
f = f + ae_cdt(c,d,t);

% Scale.
ar = R_c(c)*T*f;