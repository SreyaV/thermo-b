function [Psat rf rg] = Saturation_iT(ispecies,T)
% Return the saturation pressure Psat (Pa) and densities rf, rg (kg/m3)
% for any species i and temperature T (K).
% C.F. Edwards, 2-3-10

global toler
global Tcrit_i rcrit_i Pcrit_i 
global Ptrip_i rftrip_i rgtrip_i

% Shorten the names.
Pcrit = Pcrit_i(ispecies);
Tcrit = Tcrit_i(ispecies);
rcrit = rcrit_i(ispecies);
Ptrip  = Ptrip_i(ispecies);
rftrip = rftrip_i(ispecies);
rgtrip = rgtrip_i(ispecies);

% Set the tolerance for the convergence of the pressure.  Use the critical
% pressure as the representative value.
Ptoler = toler;
ftoler = 200; % J/kg

% Check for being "at" the critical point.
if(abs(T-Tcrit)/T < toler)
    Psat = Pcrit;
    rf = rcrit;
    rg = rcrit;
    return
end

% Check for above critical.  Return zeros for the densities in this case.
if(T > Tcrit)
    Psat = 0;
    rf = 0;
    rg = 0;
    return
end

% Find the spinodal densities and pressures at this temperature.
% If the spinodal does not exist, set the pressure to zero.
rgs = Vapor_Spinodal_iT(ispecies,T);
if(rgs == 0)
    Pgs = 0;
else
    Pgs = P_irT(ispecies,rgs,T);
end

rfs = Liquid_Spinodal_iT(ispecies,T);
if(rfs == 0)
    Pfs = 0;
else
    Pfs = P_irT(ispecies,rfs,T);
end

% Set the pressure limits within which to search at the two spinodals.
% If there is no liquid spinodal or the pressure is negative, 
% use the minimum pressure.
% If there is no vapor spinodal, use the critical pressure.
% The lack of a spinodal is signaled by the pressure being set to zero.
if((rgs == 0)||(Pgs > Pcrit))
    Pmax = Pcrit;
    rlmax = rftrip;
    rvmax = rgtrip/10;
else
    Pmax = Pgs;
    rlmax = rftrip;
    rvmax = rgs;
end

if((Pfs <= 0)||(rfs == 0)||(Pfs > Pcrit))
    Pmin = Ptrip/2;
    rlmin = rftrip;
    rvmin = rgtrip/10;
else
    Pmin = Pfs;
    rlmin = rfs;
    rvmin = rgtrip/10;
end
Pmid = (Pmin + Pmax)/2;
rlmid = (rlmin + rlmax)/2;
rvmid = (rvmin + rvmax)/2;

rl = rl_iTP(ispecies,T,Pmax,rlmax);
rv = rv_iTP(ispecies,T,Pmax,rvmax);
fmax = g_irT(ispecies,rl,T) - g_irT(ispecies,rv,T);
if(fmax == 0)
    Psat = Pmax;
    rf   = rl;
    rg   = rv;
    return
end

rl = rl_iTP(ispecies,T,Pmid,rlmid);
rv = rv_iTP(ispecies,T,Pmid,rvmid);
fmid = g_irT(ispecies,rl,T) - g_irT(ispecies,rv,T);
if(fmid == 0)
    Psat = Pmid;
    rf   = rl;
    rg   = rv;
    return
end

rl = rl_iTP(ispecies,T,Pmin,rlmin);
rv = rv_iTP(ispecies,T,Pmin,rvmin);
fmin = g_irT(ispecies,rl,T) - g_irT(ispecies,rv,T);
if(fmin == 0)
    Psat = Pmin;
    rf   = rl;
    rg   = rv;
    return
end

% Refine it by inverse quadratic interpolation.  
% See discussion of Brent's method in Numerical Recipes in C++, Cambridge
% University Press, 2002 (reprinted in 2005 with corrections).

imax = 50;
Plast = Pmid;
Pnewlast = Pmid;
for i=1:1:imax
    RR = fmid/fmax;
    SS = fmid/fmin;
    TT = fmin/fmax;
    PP = SS*(TT*(RR-TT)*(Pmax-Pmid) - (1-RR)*(Pmid-Pmin));
    QQ = (TT-1)*(RR-1)*(SS-1);
    
    if(QQ ~= 0)                         % Calculate new value if not singular.
        P = Pmid + PP/QQ;               % In bounds, new, so keep it.
        if((P >= Pmin)&&(P <= Pmax)&&(abs((P-Plast)/P) > toler))    
            Pnew = P;
        else                            % Out of bounds, so bisect.
            if(fmin*fmid > 0)
                Pnew = (Pmid + Pmax)/2;
            else
                Pnew = (Pmid + Pmin)/2;
            end
        end
        Plast = P;
    else                                % Singular, so bisect.
        if(fmin*fmid > 0)
            Pnew = (Pmid + Pmax)/2;
        else
            Pnew = (Pmid + Pmin)/2;
        end
    end
    rl = rl_iTP(ispecies,T,Pnew,rl);
    rv = rv_iTP(ispecies,T,Pnew,rv);
    gl = g_irT(ispecies,rl,T);
    gv = g_irT(ispecies,rv,T);
    
    fnew = gl - gv;
    if(fnew == 0)
        Psat = Pnew;
        rf   = rl;
        rg   = rv;
        return
    end

    f_small_enough = (abs(fnew) < ftoler);
    % Converged at successive solutions--faster
    P_close_enough = (abs((Pnew-Pnewlast)/Pnew) < Ptoler);
    % Boundaries around are converged--more secure
%     P_close_enough = (abs((Pmax-Pmin)/Pnew) < Ptoler);

    if(f_small_enough && P_close_enough)    % Good enough yet?
%         Iterations = i
        Psat = Pnew;
        rf = rl;
        rg = rv;
        break                               % Jump out!
    end

    if(Pnew > Pmid)
        Pmin = Pmid;               % Move lower bound in.
        fmin = fmid;
    else
        Pmax = Pmid;               % Move upper bound in.
        fmax = fmid;
    end
    Pmid = Pnew;
    fmid = fnew;
    Pnewlast = Pnew;
end
if(i == imax)
    Pnewlast
    Pnew
    fnew
    disp('Fell of end of interpolation loop in Saturation_iT')
    return
end
