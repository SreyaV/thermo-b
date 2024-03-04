function [T q V y x rg rf] = Flash_zhP(z,h,P,varargin)
% Return the temperature (K), quality, molar quality, compositions, and densities (kg/m3)
% for the flash problem with specified overall composition z at pressure P
% (Pa) and enthalpy h (J/kg).
% C.F. Edwards, 2-20-10

global toler

% Use information supplied if available.
switch nargin
    case 3
        % Must find a close initial point.  Get the enthalpies for bubble and
        % dew first.  Then use interpolation.
        [Tdew rgdew rfdew xdew] = Dew_cP(z,P);
        hdew = h_crT(z,rgdew,Tdew);
        [Tbub rfbub rgbub ybub] = Bubble_cP(z,P);
        hbub = h_crT(z,rfbub,Tbub);
        H = (h-hbub)/(hdew-hbub);
        % Use linear interpolation based on the enthalpy to set the starting
        % temperature for the NR procedure in the Flash_zTP calculation.
        T  = Tbub + (Tdew-Tbub)*H;
        y = zeros(1,length(z));
        for i=1:1:length(z)
            y(i)  = ybub(i) + (z(i)-ybub(i))*H;
        end
        rf = rfbub + (rfdew-rfbub)*H;
        rg = rgbub + (rgdew-rgbub)*H;
        V  = H;
    case 8
        V  = varargin{1};
        % Keep guess in bounds.
        if V >= 1
            V = 1-10*toler;
        end
        if V <= 0
            V = 0+10*toler;
        end
        T  = varargin{2};
        y  = varargin{3};
        rf = varargin{4};
        rg = varargin{5};
    otherwise
        disp('Incorrect number of arguments in Flash_zhP')
        return
end

Tinc = T*1e-5;      % Use a small increment for y derivatives.
htoler = 1;         % Use a loose absolute tolerance for h.
Ttoler = toler;     % Use a normal convergence tolerance for T.
Tlast = T;          % Save last value for convergence on y.

imax = 10;
for i=1:1:imax
    T
    % See how we are doing on T.
    [q V y x rg rf] = Flash_zTP(z,T,P,V,y,rf,rg);
    hf = h_crT(x,rf,T);
    hg = h_crT(y,rg,T);
    hest = (1-q)*hf + q*hg;
    f = hest - h;

    if((abs(f) < htoler)&&((abs(T-Tlast)/T) < Ttoler))
        % Values are set.  Go back.
        return
    end
    Tlast = T;
    Vlast = V;
    
    % Use NR to adjust temperature to get h correct.
    % Use a central difference for the numerical derivative.
    % Go low.
    Tlow = T - Tinc;
    [q V y x rg rf] = Flash_zTP(z,Tlow,P,V,y,rf,rg);
    hf = h_crT(x,rf,Tlow);
    hg = h_crT(y,rg,Tlow);
    hlow = (1-q)*hf + q*hg;
    % Go high.
    Thigh = T + Tinc;
    [q Vhigh yhigh x rghigh rfhigh] = Flash_zTP(z,Thigh,P,V,y,rf,rg);
    hf = h_crT(x,rf,Thigh);
    hg = h_crT(y,rg,Thigh);
    hhigh = (1-q)*hf + q*hg;
    
    dfdT = ((hhigh-h)-(hlow-h))/(2*Tinc);
    T = T - f/dfdT;

    % Intepolate the starting values using the new T.
    V = V + (Vhigh-V)*(T-Tlast)/(Thigh-Tlast);
    % Bisect the V guess if incorrect.
    if V > 1
        V = (1+Vlast)/2;
    end
    if V < 0
        V = (0+Vlast)/2;
    end
    for j=1:1:length(z)
        y(j)  = y(j) + (yhigh(j)-y(j))*(T-Tlast)/(Thigh-Tlast);
    end
    rg = rg + (rghigh-rg)*(T-Tlast)/(Thigh-Tlast);
    rf = rf + (rfhigh-rf)*(T-Tlast)/(Thigh-Tlast);    
end

disp('Fell off the end of NR loop for T in Flash_zhP')
