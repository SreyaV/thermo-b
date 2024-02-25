function [Tsat Psat rgsat] = Saturation_irf(i,rf)
% Return the saturation data for any given saturated liquid density.
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global Ttrip_i rftrip_i
global toler

% Check limits.
if(rf < rcrit_i(i))
    i
    disp('P out of bounds in Saturation_irf')
    return
end
if(rf > rftrip_i(i))
    i
    disp('P out of bounds in Saturation_irf')
    return
end

% Rough out a bracket that spans the solution.
dT = (Tcrit_i(i)-Ttrip_i(i))/10;
for T=Ttrip_i(i):dT:Tcrit_i(i)
    [Psat rfsat rgsat] = Saturation_iT(i,T);
    if(rfsat < rf)
        break
    end
    Rlow = rf - rfsat;
end
Tlow  = T - dT;
Thigh = T;
Rhigh = rf - rfsat;

% Use false position to hunt it down.
for n=1:1:10
    Tnext = Tlow + (Thigh - Tlow)*(-Rlow)/(Rhigh-Rlow);
    [Psat rfsat rgsat] = Saturation_iT(i,Tnext);
    Rnext = rf - rfsat;
    if(abs(Rnext/rf) <= toler)
        Tsat = Tnext;
        return
    end
    if(Rnext > 0)
        Thigh = Tnext;
        Rhigh = Rnext;
    else
        Tlow  = Tnext;
        Rlow  = Rnext;
    end
end
disp('Slow convergence by False Position in Saturation_irf.  Try bisection...')

% If you still don't have it, try bisection.
for n=1:1:100
    Tmid = (Tlow + Thigh)/2;
    [Psat rfsat rgsat] = Saturation_iT(i,Tmid);
    Rmid = rf - rfsat;
    if(abs((Thigh - Tlow)/Tmid) <= toler)
        Tsat = Tmid;
        return
    end
    if(Rmid > 0)
        Thigh = Tmid;
    else
        Tlow  = Tmid;
    end
end    
disp('Fell of end of loop in Saturation_irf')
