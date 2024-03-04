function [Tcrit rcrit] = Critical_i(ispecies)
% Find the numerical critical temperature (K) and density (kg/m3) for a
% given species i.
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global toler

% Use 2x 1-D Newton-Raphson to find the T-rho location where both
% derivatives are zero.  (A plot of the contours of the two derivatives
% near the experimental critical point will show why you can do this.)
imax = 20;
jmax = 20;
kmax = 50;
Tred = Tcrit_i(ispecies);
rred = rcrit_i(ispecies);
T = Tred;
r = rred;
for k=1:1:kmax
    % Do N-R on the first derivative to find a T that gives zero.
    Tlast = T;
    for i=1:1:imax
        f = dPdr_irT(ispecies,r,T);
        dfdT = d2PdrdT_irT(ispecies,r,T);
        T = T - f/dfdT;
        if((f < toler/10)&&(abs(T-Tlast)/T < toler/10))
            break
        end
        Tlast = T;
    end
    % Do N-R on the second derivative to find a rho that gives zero.
    rlast = r;
    for j=1:1:jmax
        g = d2Pdr2_irT(ispecies,r,T);
        dgdr = d3Pdr3_irT(ispecies,r,T);
        r = r - g/dgdr;
        if((g < toler/10)&&(abs(r-rlast)/r < toler/10))
            break
        end
        rlast = r;
    end
    % Both criteria must be satisfied to exit.
    f = dPdr_irT(ispecies,r,T);
    if(abs(f) < toler/10)
        break
    end
end

% Check to see if you fell through the loops.
if(k == kmax)
    ispecies
    r
    T
    disp('Fell off end of NR loop for f & g in Critical_Point_i')
end
Tcrit = T;
rcrit = r;
