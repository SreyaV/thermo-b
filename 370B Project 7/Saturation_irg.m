function [Tsat Psat rfsat] = Saturation_irg(i,rg)
% Return the saturation data for any given saturated vapor density.
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i
global Ttrip_i rgtrip_i
global toler

% Check limits.
if(rg > rcrit_i(i))
    i
    disp('P out of bounds in Saturation_irg')
    return
end
if(rg < rgtrip_i(i))
    i
    disp('P out of bounds in Saturation_irg')
    return
end

% Rough out a bracket that spans the solution.
dT = (Tcrit_i(i)-Ttrip_i(i))/10;
for T=Ttrip_i(i):dT:Tcrit_i(i)
    [Psat rfsat rgsat] = Saturation_iT(i,T);
    if(rgsat > rg)
        break
    end
    Rlow = rgsat - rg;
end
Tlow  = T - dT;
Thigh = T;
Rhigh = rgsat - rg;

% Use false position to hunt it down.
for n=1:1:10
    Tnext = Tlow + (Thigh - Tlow)*(-Rlow)/(Rhigh-Rlow);
    [Psat rfsat rgsat] = Saturation_iT(i,Tnext);
    Rnext = rgsat - rg;
    if(abs(Rnext/rg) <= toler)
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
disp('Slow convergence by False Position in Saturation_irg.  Try bisection...')

% If you still don't have it, try bisection.
for n=1:1:100
    Tmid = (Tlow + Thigh)/2;
    [Psat rfsat rgsat] = Saturation_iT(i,Tmid);
    Rmid = rgsat - rg;
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
disp('Fell of end of loop in Saturation_irg')
