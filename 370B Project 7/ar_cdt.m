function ar = ar_cdt(c,d,t)
% Return the normalized residual Helmholtz function for composition c,
% reduced density d, and reduced temperature t.

% Make an ideal solution.
f = 0;
for(i=1:1:length(c))
    if(c(i) ~= 0)
        f = f + c(i)*(ar_idt(i,d,t) + log(c(i)));
    end
end

% Add the excess.
f = f + ae_cdt(c,d,t);

% Scale.
ar = f;