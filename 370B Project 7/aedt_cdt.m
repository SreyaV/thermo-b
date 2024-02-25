function aedt = aedt_cdt(c,d,t)
% Return the mixed derivative of the excess contribution wrt delta and tau
% at a specified value of composition (c), delta (d), and tau (t).

global Fmix Mix_i Mix_j Mix_N NMix

% Use the expression for the fit.

f = 0;
for(k=1:1:NMix)
    f = f + Mix_N(k)*Mix_i(k)*d^(Mix_i(k)-1)*Mix_j(k)*t^(Mix_j(k)-1);
end

aedt = 0;
for(i=1:1:length(c)-1)
    for(j=i+1:1:length(c))
        aedt = aedt + c(i)*c(j)*Fmix(i,j)*f;
    end
end