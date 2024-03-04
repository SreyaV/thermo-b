function [Tinfl rinfl] = Pr_Inflection_c(c)
% Calculate the inflection point temperature (K) and density (kg/m3) for a
% given composition (molefractions).
% C.F. Edwards, 1/15/09

global toler

% Check to make sure you have a viable composition.
csum = sum(c);
if(abs(csum-1) > toler)
    c
    disp('Molefractions do not sum to unity in Pr_Inflection')
    return
end
if(~isreal(c))
    c
    disp('Molefractions have imaginary components in Pr_Inflection')
    return
end
if(max(c) > 1)
    c
    disp('Component molefraction exceeds unity in Pr_Inflection')
    return
end
if(min(c) < 0)
    c
    disp('Component molefraction less than zero in Pr_Inflection')
    return
end

% Remove any trace of bad normalization.
c = c/csum;

% Check for pure components.
for i=1:1:length(c)
    if(abs(c(i)-1) < toler)
        % Do you want the experimental values?
%         Tinfl = Tcrit_i(i);
%         rinfl = rcrit_i(i);
        % Or do you want the actual inflection point?
        [Tinfl rinfl] = Critical_i(i);
        return
    end
end

% Get the reducing values.  These are used to provide a starting point.
Tred = Tred_c(c);
rred = rred_c(c);

% Solve by using two 1-D N-R sweeps in succession.
imax = 20;
jmax = 20;
kmax = 50;
T = Tred;
r = 0.9*rred;
for k=1:1:kmax
    Tlast = T;
    % N-R of f = dPdr in the T direction.
    for i=1:1:imax
        f = dPdr_crT(c,r,T);
        dfdT = d2PdrdT_crT(c,r,T);
        % Limit the step size.
        mag_dT = abs(f/dfdT);
        sign_dT = (f/dfdT)/mag_dT;
        if(mag_dT > 1)
            dT = sign_dT*1;
        else
            dT = f/dfdT;
        end
        T = T - dT;
        % Close enough?  Check magnitude and bound.
        if((abs(f) < toler)&&(abs(T-Tlast)/T < toler))
            break
        end
        Tlast = T;
    end
    i_Iterations = i;
    
    % N-R of g = d2Pdr2 in the r direction.
    rlast = r;
    for j=1:1:jmax
        g = d2Pdr2_crT(c,r,T);
        dgdr = d3Pdr3_crT(c,r,T);
        % Limit the step size.
        mag_dr = abs(g/dgdr);
        sign_dr = (g/dgdr)/mag_dr;
        if(mag_dr > 10)
            dr = sign_dr*1;
        else
            dr = g/dgdr;
        end
        r = r - dr;
        % Close enough?  Check magnitude and bound.
        if((abs(g) < toler)&&(abs(r-rlast)/r < toler))
            break
        end
        rlast = r;
    end
    j_Iterations = j;
    
    % We know g is converged.  Is f still small enough?
    f = dPdr_crT(c,r,T);
    if(abs(f) < toler)
        Tinfl = T;
        rinfl = r;
        return
    end
end
disp('Fell off end of NR loop for f & g in Pr_Inflection_c')
