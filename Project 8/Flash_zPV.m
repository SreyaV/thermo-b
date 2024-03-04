function [T y x rg rf] = Flash_zPV(z,P,V,Tstart,ystart,rlstart,rvstart,direction)
% Return the temperature, compositions, and densities (kg/m3) for the flash
% problem with specified overall composition z at pressure P (Pa) and with
% molar quality V (mole fraction vapor).

global toler
global N2 O2 Ar

if(direction ~= 0)
    T  = Tstart;
    y  = ystart;
    rf = rlstart;
    rg = rvstart;
else
    [Tdew rgdew rfdew xdew] = Dew_Point_cP(z,P,0,0,0,0,0);
    Tdew;
    [Tbub rfbub rgbub ybub] = Bubble_Point_cP(z,P,0,0,0,0,0);
    Tbub;
    % Use linear interpolation based on the quality to set the starting
    % temperature for the NR procedure in the Flash_zTP calculation.
    T  = Tbub + (Tdew-Tbub)*V;
    for(i=1:1:length(z))
        y(i)  = ybub(i) + (z(i)-ybub(i))*V;
    end
    rf = rfbub + (rfdew-rfbub)*V;
    rg = rgbub + (rgdew-rgbub)*V;
end

Tinc = T*1e-5;      % Use a small increment for y derivatives.
Vtoler = 10*toler;  % Use a looser tolerance for V than T.
Ttoler = toler;     % Use normal tolerance for T.
Tlast = T;          % Save last value for convergence on y.

imax = 10;
for(i=1:1:imax)
    T;
    % See how we are doing on T.
    [q Vflash y x rg rf] = Flash_zTP(z,T,P,V,y,rf,rg,1);
    Vflash;
    f = Vflash - V;

    if((abs(f/V) < Vtoler)&&((abs(T-Tlast)/T) < Ttoler))
        break
    end
    Tlast = T;
    
    % Use NR to adjust temperature to get V correct.
    % Use a central difference for the numerical derivative.
    Tlow = T - Tinc;
    [q Vlow y x rg rf] = Flash_zTP(z,Tlow,P,V,y,rf,rg,1);
    Vlow;
    
    Thigh = T + Tinc;
    [q Vhigh y x rg rf] = Flash_zTP(z,Thigh,P,V,y,rf,rg,1);
    Vhigh;
    
    dfdT = ((Vhigh-V)-(Vlow-V))/(2*Tinc);
    T = T - f/dfdT;
end
iT = i;
if(i == imax)
    disp('Fell off the end of NR loop for T in Flash_zPV')
end
T = T;
y = y;
x = x;
rg = rg;
rf = rf;