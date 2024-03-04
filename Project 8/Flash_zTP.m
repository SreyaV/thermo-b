function [q V y x rg rf] = Flash_zTP(z,T,P,varargin)
% Return the qualities, compositions, and densities (kg/m3) for the flash
% problem with specified overall composition z at temperature T (K) and 
% pressure P (Pa).
% C.F. Edwards, 2-20-10

Ptoler = 1.5e-4;

% Use information supplied if available.
switch nargin
    case 3
        % Get the bubble and dew points as quality limits.
        [Pdew rgdew rfdew xdew] = Dew_cT(z,T);
        [Pbub rfbub rgbub ybub] = Bubble_cT(z,T);
        % Check for the low end.
        if(abs(P-Pbub)/Pbub < Ptoler)
            q = 0;
            V = 0;
            y = ybub;
            x = z;
            rg = rgbub;
            rf = rfbub;
            return
        end
        % Check for the high end.
        if(abs(P-Pdew)/Pdew < Ptoler)
            q = 1;
            V = 1;
            y = z;
            x = xdew;
            rg = rgdew;
            rf = rfdew;
            return
        end
        % Check that we are in between.
        if((P > Pbub)||(P < Pdew))
            disp('P out of bounds in Flash_zTP')
            return
        end
        % Use linear interpolation based on the pressure to set the starting
        % point for the NR procedure for V.
        V = 1 - (P-Pdew)/(Pbub-Pdew);
        y = zeros(1,length(z));
        for i=1:1:length(z)
            y(i) = ybub(i) + V*(z(i)-ybub(i));
        end
        rl = rfbub;
        rv = rgdew;
    case 7
        V  = varargin{1};
        y  = varargin{2};
        rl = varargin{3};
        rv = varargin{4};
    otherwise
        disp('Incorrect number of arguments in Flash_zTP')
        return
end

ftoler = 1;     % Tolerance for convergence of chemical potential on isp.
Vtoler = 1e-5;  % Tolerance for convergence of molar quality.
Vinc   = 1e-5;  % Use a small increment for the molar quality derivatives.
Vlast  = V;     % Save last value for convergence on V.

% Find the largest species.  Use for closure.
[maxsp isp] = max(z);

imax = 200;
for i=1:1:imax
    % Find the matching liquid compositions for all but isp.
    [y x] = yx_izTPV(isp,z,T,P,V,y,rl,rv);
    
    % See how we are doing on isp.
    rv = rv_cTP(y,T,P,rv);
    muiv = mui_icrT(isp,y,rv,T);
    rl = rl_cTP(x,T,P,rl);
    muil = mui_icrT(isp,x,rl,T);
    fisp = (muiv - muil);

    if((abs(fisp) < ftoler)&&((abs(V-Vlast)/V) < Vtoler))
        rg = rv;
        rf = rl;
        % Get the mass-based quality.
        Mv = M_c(y);
        Ml = M_c(x);
        q = Mv*V/(Mv*V + Ml*(1-V));
        return
    end
    Vlast = V;
    
    % Use NR to adjust the molar quality to get isp correct.

    % Use a central difference for the numerical derivative.
    % Go low.
    [ylow xlow] = yx_izTPV(isp,z,T,P,V-Vinc,y,rl,rv);
    rvlow = rv_cTP(ylow,T,P,rv);
    muivlow = mui_icrT(isp,ylow,rvlow,T);
    rllow = rl_cTP(xlow,T,P,rl);
    muillow = mui_icrT(isp,xlow,rllow,T);
    % Go high.
    [yhigh xhigh] = yx_izTPV(isp,z,T,P,V+Vinc,y,rl,rv);
    rvhigh = rv_cTP(yhigh,T,P,rv);
    muivhigh = mui_icrT(isp,yhigh,rvhigh,T);
    rlhigh = rl_cTP(xhigh,T,P,rl);
    muilhigh = mui_icrT(isp,xhigh,rlhigh,T);
    
    dfispdV = ((muivhigh-muilhigh)-(muivlow-muillow))/(2*Vinc);
    dV = -fisp/dfispdV;
    % Limit the step size.  This can help if N-R tries to run off.
%     limit = 0.05;
%     if(abs(dV) > limit)
%         dV = sign(dV)*limit;
%     end
    V = V + dV;
    
    % Check the limits.  Bisect if out of bounds.
    if(V < 0)
        V = (Vlast-0)/2;
    end
    if(V > 1)
        V = (1-Vlast)/2;
    end
end

% Allow the procedure to fall off the end and still return a value.
disp('Fell off the end of NR loop for y in Flash_zTP')
format long
V
Vlast
format short
rg = rv;
rf = rl;
% Get the mass-based quality.
Mv = M_c(y);
Ml = M_c(x);
q = Mv*V/(Mv*V + Ml*(1-V));
