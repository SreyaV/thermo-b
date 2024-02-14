function [r T rf rg] = rT_ihP(ispecies,h,P)
% Return the density (kg/m3) and T (K) for given h (J/kg) and P (Pa).
% If the state point is under the dome, return the saturation densities
% too.  If not, return zeros in these variables.
% C.F. Edwards, 2-6-11

global Tcrit_i rcrit_i Ttrip_i Tupper_i
global toler
persistent inumcrit Tnumcrit rnumcrit Pnumcrit

% The first time through you must set the species.
if(isempty(inumcrit))
    inumcrit = ispecies;
end

% Get numerical critical point if not already available.
if(isempty(Pnumcrit)||(inumcrit ~= ispecies))
    inumcrit = ispecies;
    try
        % This will work if you have Critical_i
        [Tnumcrit rnumcrit] = Critical_i(ispecies);
    catch
        % If you don't it is Project 6 on nH2.
        disp('Species must be nH2 in rT_ihP')
        Tnumcrit = 32.93800;
        rnumcrit = 31.35997;
    end
    Pnumcrit = P_irT(inumcrit,rnumcrit,Tnumcrit);
end

% Shorten the names.
Tcrit  = Tcrit_i(ispecies);
rcrit  = rcrit_i(ispecies);
Ttrip  = Ttrip_i(ispecies);
Tupper = Tupper_i(ispecies);
hcrit  = h_irT(ispecies,rcrit,Tcrit);

% Walk along an isobar incrementing temperature to find the enthalpy.
% Watch out for the saturation region.

Nints = 10;                 % Number of roughing intervals.

if(P > Pnumcrit)
    % The pressure is supercritical.  No need to worry about the dome.
    Thigh = Tupper;         % Set search bounds.
    Tlow  = Ttrip;
    vapor = 0;              % Logic flag for density
elseif(h > 2*hcrit)
    % The enthalpy is high enough to not worry about the dome.
    Thigh = Tupper;         % Set search bounds.
    Tlow  = Tcrit;
    vapor = 1;              % Logic flag for density    
else
    % The pressure is below critical so we need to look at saturation.
    [Tsat rf rg] = Saturation_iP(ispecies,P);
    hf   = h_irT(ispecies,rf,Tsat);
    hg   = h_irT(ispecies,rg,Tsat);
    
    % Branch depending on where the enthalpy lies.
    if((h >= hf)&&(h <= hg))
        % We are under the dome.  Use quality.
        q = (h - hf)/(hg - hf);
        v = 1/rf + q*(1/rg - 1/rf);
        r = 1/v;
        T = Tsat;
        return
    end
    
    % Not under dome.  Either above or below.
    if(h > hg)
        Thigh = Tupper;     % Superheated vapor.
        Tlow  = Tsat;
        vapor = 1;          % Logic flag for density
    end
    if(h < hf)
        Thigh = Tsat;       % Subcooled liquid.
        Tlow  = Ttrip;
        vapor = 0;          % Logic flag for density
    end
end

% The bounds are set and the saturation region is elminated.
% Walk the temperature from low to high at constant pressure.

dT = (Thigh-Tlow)/Nints;
% Watch out for out of bounds on low T end.
hmin = h_irT(ispecies,rl_iTP(ispecies,Ttrip,P),Ttrip);
hlast = hmin;
for T=Tlow:dT:Thigh
    if(vapor)
        rnext = rv_iTP(ispecies,T,P);
    else
        rnext = rl_iTP(ispecies,T,P);
    end
    hnext = h_irT(ispecies,rnext,T);
    if(hnext >= h)
        break       % Crossed target enthalpy.
    else
        hlast = hnext;
    end
end
Thigh = T;
Tlow  = T-dT;
Reshigh = hnext - h;
Reslow  = hlast - h;
if(Reslow*Reshigh > 0)
    disp('Failed to bracket in rT_hP.')
    if hlast==hmin
        disp('Enthapy requires temperature below triple point...')
        disp(' ')
    end
    return
end

% Start a false position loop to run this to tolerance.

NFPs = 100;
for i=1:1:NFPs
    Tnext = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow);
    if(vapor)
        rnext = rv_iTP(ispecies,Tnext,P);
    else
        rnext = rl_iTP(ispecies,Tnext,P);
    end
    hnext = h_irT(ispecies,rnext,Tnext);
    Res   = hnext - h;
    if(abs(Res/h) < toler)
        T = Tnext;          % Found the solution to tolerance.
        r = rnext;
        rf = 0;
        rg = 0;
        iterations = i;
        return              % Send r,T back.
    end
    if(Res > 0)
        Thigh   = Tnext;    % Midpoint data is on the high side of P.
        Reshigh = Res;
    else
        Tlow    = Tnext;    % Midpoint data is on the low side of P.
        Reslow  = Res;
    end
end
disp('Failed to converge in FP loop in rT_ihP')
vapor
Tnext
hnext
i