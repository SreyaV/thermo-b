function [T x y rf rg] = Flash_zqP(z,q,P,varargin)
% Return the temperature, compositions, and densities (kg/m3) for the flash
% problem with specified overall composition z at pressure P (Pa) and with
% quality q (mass fraction vapor).
% C.F. Edwards, 2-20-10

global toler

% Use information supplied if available.
switch nargin
    case 3
        % Get the bubble and dew points (q = 0 and 1 respectively).
        [Tbub rfbub rgbub ybub] = Bubble_cP(z,P);
        [Tdew rgdew rfdew xdew] = Dew_cP(z,P);
        % Use linear interpolation based on the quality to set the starting
        % temperature for the NR procedure in the Flash_zTP calculation.
        T  = Tbub + (Tdew-Tbub)*q;
        y = zeros(1,length(z));
        for i=1:1:length(z)
            y(i)  = ybub(i) + (z(i)-ybub(i))*q;
        end
        V = q*M_c(z)/M_c(y);            % Molar quality
        rf = rfbub + (rfdew-rfbub)*q;
        rg = rgbub + (rgdew-rgbub)*q;
    case 7
        T  = varargin{1};
        y  = varargin{2};
        rf = varargin{3};
        rg = varargin{4};
        V  = q*M_c(z)/M_c(y);           % Molar quality
    otherwise
        disp('Incorrect number of arguments in Flash_zqP')
        return
end

Tinc = T*1e-6;      % Use a small increment for y derivatives.
qtoler = 1e-4;      % Use a looser tolerance for q than T.
Ttoler = toler;     % Use normal tolerance for T.
Tlast = T;          % Save last value for convergence on y.

imax = 50;
for i=1:1:imax
    % See how we are doing on T.
    if(P > 5e4)
        [qflash V y x rg rf] = Flash_zTP(z,T,P,V,y,rf,rg);
    else
        [qflash V y x rg rf] = Flash_zTP(z,T,P);
    end
    f = qflash - q;

    if((abs(f) < qtoler)&&((abs(T-Tlast)/T) < Ttoler))
        break
    end
    Tlast = T;
    Vlast = V;
    
    % Use NR to adjust temperature to get q correct.
    % Use a central difference for the numerical derivative.
    [qlow V y x rg rf] = Flash_zTP(z,T-Tinc,P,V,y,rf,rg);
    [qhigh V y x rg rf] = Flash_zTP(z,T+Tinc,P,V,y,rf,rg);
    dfdT = ((qhigh-q)-(qlow-q))/(2*Tinc);

%     % Use a forward difference for the numerical derivative.
%     [qhigh V y x rg rf] = Flash_zTP(z,T+Tinc,P,V,y,rf,rg);
%     dfdT = (qhigh-qflash)/Tinc;

    % Update the estimate of T.
    dT = -f/dfdT;
    T = T + dT;
    V = Vlast;
end

if(i == imax)
    V
    disp('Fell off the end of NR loop for T in Flash_zqP')
end

