function [r T q V y x rg rf] = rT_chP(c,h,P)
% Return the density (kg/m3) and T (K) for given h (J/kg) and P (Pa).
% If the state point is under the dome, return the tie-line information
% too.  If not, return zeros in these variables.

global toler
persistent cinfl Tinfl rinfl Pinfl

% The first time through you must set the composition.
% Get numerical critical point if not already available.
if(isempty(cinfl))
    cinfl = c;
    [Tinfl rinfl] = Pr_Inflection_c(c);
    Pinfl = P_crT(cinfl,rinfl,Tinfl);
elseif(cinfl ~= c)
    cinfl = c;
    [Tinfl rinfl] = Pr_Inflection_c(c);
    Pinfl = P_crT(cinfl,rinfl,Tinfl);    
end

% Set some boundaries.
Tupper = 870;
Tlower = 60;
hinfl  = h_crT(cinfl,rinfl,Tinfl);

% Walk along an isobar incrementing temperature to find the enthalpy.
% Watch out for the saturation region.

Nints = 10;                 % Number of roughing intervals.

if(P > Pinfl)
    % The pressure is above inflection.  No need to worry about the dome.
    Thigh = Tupper;         % Set search bounds.
    Tlow  = Tlower;
    vapor = 0;              % Logic flag for density
elseif(h > 100e3)
    % The enthalpy is high enough to not worry about the dome.
    Thigh = Tupper;         % Set search bounds.
    Tlow  = Tinfl;
    vapor = 1;              % Logic flag for density    
else
    % The pressure is below inflection so we need to look at saturation.
    [Tfsat rfsat rg y] = Bubble_cP(c,P);
    [Tgsat rgsat rf x] = Dew_cP(c,P);
    hf   = h_crT(c,rfsat,Tfsat);
    hg   = h_crT(c,rgsat,Tgsat);
    
    % Check for sitting right at the edge of the dome.  This can happen in
    % distillation problems.
    if(abs(h-hg)/(hg-hf) < 1e-4)
        % We are at saturated vapor.
        T  = Tgsat;
        q  = 1;
        V  = 1;
        y  = c;
        rg = rgsat;
        r  = rg;
        return
    elseif(abs(h-hf)/(hg-hf) < 1e-4)
        % We are at saturated liquid.
        T  = Tfsat;
        q  = 0;
        V  = 0;
        x  = c;
        rf = rfsat;
        r  = rf;
        return
    end
    
    % Branch depending on where the enthalpy lies.
    if((h >= hf)&&(h <= hg))
        % We are under the dome.  Flash it.
        [T q V y x rg rf] = Flash_zhP(c,h,P);
        v = 1/rf + q*(1/rg - 1/rf);
        r = 1/v;
        return
    end
    
    % Not under dome.  Either above or below.
    if(h > hg)
        Thigh = Tupper;     % Superheated vapor.
        Tlow  = Tgsat;
        vapor = 1;          % Logic flag for density
    end
    if(h < hf)
        Thigh = Tfsat;       % Subcooled liquid.
        Tlow  = Tlower;
        vapor = 0;          % Logic flag for density
    end
end

% The bounds are set and the saturation region is elminated.
% Walk the temperature from low to high at constant pressure.

dT = (Thigh-Tlow)/Nints;
for T=Tlow:dT:Thigh
    if(vapor)
        rnext = rv_cTP(c,T,P);
    else
        rnext = rl_cTP(c,T,P);
    end
    hnext = h_crT(c,rnext,T);
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
    disp('Failed to bracket in rT_chP')
    return
end

% Start a false position loop to run this to tolerance.

NFPs = 100;
for i=1:1:NFPs
    Tnext = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow);
    if(vapor)
        rnext = rv_cTP(c,Tnext,P);
    else
        rnext = rl_cTP(c,Tnext,P);
    end
    hnext = h_crT(c,rnext,Tnext);
    Res   = hnext - h;
    if(abs(Res/h) < toler)
        T = Tnext;          % Found the solution to tolerance.
        r = rnext;
        q = 0;
        V = 0;
        y = zeros(length(c));
        x = y; 
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
disp('Failed to converge in FP loop in rT_chP')
vapor
Tnext
hnext
i