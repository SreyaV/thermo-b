function [Tsat rf rg] = Saturation_iP(i,P)
% Return the saturation temperature (K) data for any given saturation pressure.
% C.F. Edwards, 2-3-10

global Tcrit_i Pcrit_i
global Ttrip_i Ptrip_i
global toler

% Check limits.
if(P > Pcrit_i(i))
    i
    disp('P out of bounds in Saturation_iP')
    return
end
if(P < Ptrip_i(i))
    i
    disp('P out of bounds in Saturation_iP')
    return
end

% Rough out a bracket that spans the solution.
dT = (Tcrit_i(i)-Ttrip_i(i))/10;
for T=Ttrip_i(i):dT:Tcrit_i(i)
    [Psat rf rg] = Saturation_iT(i,T);
    if(Psat > P)
        break
    end
    Rlow = Psat - P;
end
Tlow  = T - dT;
Thigh = T;
Rhigh = Psat - P;

% Use false position to hunt it down.
for n=1:1:10
    Tnext = Tlow + (Thigh - Tlow)*(-Rlow)/(Rhigh-Rlow);
    [Psat rf rg] = Saturation_iT(i,Tnext);
    Rnext = Psat - P;
    if(abs(Rnext/P) <= toler)
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
disp('Slow convergence by False Position in Saturation_iP.  Try bisection...')

% If you still don't have it, try bisection.
for n=1:1:100
    Tmid = (Tlow + Thigh)/2;
    [Psat rf rg] = Saturation_iT(i,Tmid);
    Rmid = Psat - P;
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
disp('Fell of end of loop in Saturation_iP')
