function aett = aett_cdt(c,d,t)
% Return the second derivative of the excess contribution wrt tau for
% a specified value of composition (c), delta (d), and tau (t).

global Fmix Mix_i Mix_j Mix_N NMix

% Use the expression for the fit.

f = 0;
for(k=1:1:NMix)
    f = f + (Mix_j(k)-1)*Mix_N(k)*d^Mix_i(k)*Mix_j(k)*t^(Mix_j(k)-2);
end

aett = 0;
for(i=1:1:length(c)-1)
    for(j=i+1:1:length(c))
        aett = aett + c(i)*c(j)*Fmix(i,j)*f;
    end
end