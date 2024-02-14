function ardd = ardd_idt(i,d,t)
% Return the second derivative of the residual normalized Helmholtz
% function of species i with respect to delta at a specified value 
% of delta (d), and tau (t).
% C.F. Edwards, 2-3-10

global FR_Npoly FR_Nexp FR_Ngaus FR_Nnonan
global FR_N FR_d FR_t FR_c
global FR_eta FR_beta FR_gamma FR_epsilon
global FR_a FR_b FR_NAbeta FR_A FR_B FR_C FR_D

% Build up the residual by the types of terms.
ardd = 0;
k = 1;

% Start with the polynomial terms:
for j=1:1:FR_Npoly(i)
    ardd = ardd + (FR_d(k,i)-1)*FR_d(k,i)*FR_N(k,i)*d^(FR_d(k,i)-2)*t^FR_t(k,i);
    k = k+1;
end

% Add the exponential terms:
for j=1:1:FR_Nexp(i)
    ardd = ardd + ...
         ( FR_d(k,i)*(FR_d(k,i)-1)*d^(FR_d(k,i)-2) - ...
           FR_d(k,i)*FR_c(k,i)*d^(FR_c(k,i)-1)*d^(FR_d(k,i)-1) - ...
           d^(FR_c(k,i)+FR_d(k,i)-2)*FR_c(k,i)*(FR_c(k,i)-1) - ...
           d^(FR_c(k,i)+FR_d(k,i)-2)*FR_c(k,i)*FR_d(k,i) + ...
           d^(FR_c(k,i)+FR_d(k,i)-2)*FR_c(k,i)*FR_c(k,i)*d^FR_c(k,i) ) ... 
         * FR_N(k,i)*t^FR_t(k,i)*exp(-d^FR_c(k,i));
    k = k+1;
end

% Add the Gaussian terms:
for j=1:1:FR_Ngaus(i)
    ardd = ardd + ...
         ( 1*FR_d(k,i)*(FR_d(k,i)-1)*d^(FR_d(k,i)-2) - ...
           2*FR_d(k,i)*FR_eta(k,i)*(d-FR_epsilon(k,i))*d^(FR_d(k,i)-1) - ...
           2*FR_eta(k,i)*d^FR_d(k,i) - ...
           2*FR_d(k,i)*FR_eta(k,i)*(d-FR_epsilon(k,i))*d^(FR_d(k,i)-1) + ...
           4*FR_eta(k,i)*(d-FR_epsilon(k,i))^2*FR_eta(k,i)*d^FR_d(k,i) )...
         * FR_N(k,i)*t^FR_t(k,i)*exp(-FR_eta(k,i)*(d-FR_epsilon(k,i))^2-FR_beta(k,i)*(t-FR_gamma(k,i))^2);
    k = k+1;
end

% Add the nonanalytic terms:
for j=1:1:FR_Nnonan(i)
    Theta = (1-t) + FR_A(k,i)*((d-1)^2)^(1/(2*FR_NAbeta(k,i)));
    Delta = Theta^2 + FR_B(k,i)*((d-1)^2)^FR_a(k,i);
    Psi   = exp(-FR_C(k,i)*(d-1)^2-FR_D(k,i)*(t-1)^2);

    % Don't compute the non-analytic terms if sitting on critical.
    % Their limit as you approach critical is zero.
    if((Theta == 0)&&(Delta == 0)&&(Psi == 1))
        ardd = ardd + 0;
    else
        dTheta  = (d-1)*(1/(FR_NAbeta(k,i)))*FR_A(k,i)*((d-1)^2)^(1/(2*FR_NAbeta(k,i))-1);
        dDelta  = 2*Theta*dTheta + 2*(d-1)*FR_a(k,i)*FR_B(k,i)*((d-1)^2)^(FR_a(k,i)-1);
        dDeltab = FR_b(k,i)*Delta^(FR_b(k,i)-1)*dDelta;
        dPsi    = -2*FR_C(k,i)*(d-1)*Psi;

        d2Theta  = FR_NAbeta(k,i)^(-1)*FR_A(k,i)*((d-1)^2)^(1/(2*FR_NAbeta(k,i))-1)...
            + 2*(1/(2*FR_NAbeta(k,i))-1)*(d-1)^2*FR_NAbeta(k,i)^(-1)*FR_A(k,i)*((d-1)^2)^(1/(2*FR_NAbeta(k,i))-2);
        d2Delta  = 2*dTheta^2 + 2*Theta*d2Theta...
            + 2*FR_a(k,i)*FR_B(k,i)*((d-1)^2)^(FR_a(k,i)-1)...
            + 4*(FR_a(k,i)-1)*(d-1)^2*FR_a(k,i)*FR_B(k,i)*((d-1)^2)^(FR_a(k,i)-2);
        d2Deltab = FR_b(k,i)*(Delta^(FR_b(k,i)-1)*d2Delta + (FR_b(k,i)-1)*Delta^(FR_b(k,i)-2)*dDelta^2);
        d2Psi    = (2*FR_C(k,i)*(d-1)^2 - 1)*2*FR_C(k,i)*Psi;

        ardd = ardd + FR_N(k,i)*(Delta^FR_b(k,i)*(2*dPsi + d*d2Psi) + 2*dDeltab*(Psi + d*dPsi) + d2Deltab*d*Psi);
    end
    k = k+1;
end
