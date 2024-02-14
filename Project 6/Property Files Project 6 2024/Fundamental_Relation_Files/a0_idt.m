function a0 = a0_idt(i,d,t)
% Return the ideal gas contribution to the normalized Helmholtz function of
% species i at a specified value of delta (d), and tau (t).
% C.F. Edwards, 2-3-10

global O2
global FR_Npoly0 FR_Neinst FR_Nspec0
global FR_N0 FR_t0 FR_gamma0 FR_spec0

% Build up the residual by the types of terms.
a0 = log(d);
k = 1;

% Start with the polynomial terms:
for j=1:1:FR_Npoly0(i)
    a0 = a0 + FR_N0(k,i)*t^FR_t0(k,i);
    k = k+1;
end

% Add the logarithmic terms:
a0 = a0 + FR_N0(k,i)*log(t);
k = k+1;
a0 = a0 + FR_N0(k,i)*t*log(t);
k = k+1;

% Add the Einstein terms:
for j=1:1:FR_Neinst(i)
    a0 = a0 + FR_N0(k,i)*log(1-exp(-FR_gamma0(k,i)*t));
    k = k+1;
end

% If there are no special terms, return the answer now.
if(FR_Nspec0(i) == 0)
    return
end

% Add the special terms:
switch(i)
    case(O2)
        a0 = a0 + FR_N0(k,i)*log(exp(FR_spec0(k,i)*t)-1);
        k = k+1;
        a0 = a0 + FR_N0(k,i)*log(1+(2/3)*exp(-FR_spec0(k,i)*t));
        return
    otherwise
        disp('Mistaken attempt to execute special terms in a0_idt')
        return
end
