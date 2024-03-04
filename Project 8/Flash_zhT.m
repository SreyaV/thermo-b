function [P q V y x rg rf] = Flash_zhT(z,h,T,Vstart,Pstart,ystart,rlstart,rvstart,direction)
% Return the pressure (Pa), quality, molar quality, compositions, and densities (kg/m3)
% for the flash problem with specified overall composition z at temperature
% T (K) and enthalpy h (J/kg).

global toler
global N2 O2 Ar N

if(direction ~= 0)
    P  = Pstart;
    y  = ystart;
    V  = Vstart;
    rf = rlstart;
    rg = rvstart;
else
    % Find the bubble and dew pressures at this T.
    [Pbub rfbub rgbub ybub] = Bubble_Point_cT(z,T,0,0,0,0,0);
    [Pdew rgdew rfdew xdew] = Dew_Point_cT(z,T,0,0,0,0,0);
    
    % Divide the range between these into fsteps;
    fsteps = 20;
    dP = (Pbub-Pdew)/fsteps;
    i = 1;
    qi(i)   = 0;
    Vi(i)   = 0;
    Pi(i)   = Pbub;
    xi(i,:) = z;
    yi(i,:) = ybub;
    rfi(i)  = rfbub;
    rgi(i)  = rgbub;
    hi(i)   = h_crT(z,rfbub,T);
    V = 0;
    y = ybub;
    rf = rfbub;
    rg = rgbub;
    % Don't take the last step.  q > .9 is not stable.
    for(i=2:1:fsteps-1)
        P = Pbub - (i-1)*dP;
        [q V y x rg rf] = Flash_zTP(z,T,P,V,y,rf,rg,1);
        q
        qi(i) = q;
        Vi(i) = V;
        Pi(i) = P;
        xi(i,:) = x;
        yi(i,:) = y;
        rfi(i)  = rf;
        rgi(i)  = rg;
        hf = h_crT(x,rf,T);
        hg = h_crT(y,rg,T);
        hi(i) = (1-q)*hf + q*hg;
    end
    i = fsteps;
    qi(i)   = 1;
    Vi(i)   = 1;
    Pi(i)   = Pdew;
    xi(i,:) = xdew;
    yi(i,:) = z;
    rfi(i)  = rfdew;
    rgi(i)  = rgdew;
    hi(i)   = h_crT(z,rgdew,T);

    % Interpolate the known values to estimate the T, etc. for this h.
    P = interp1(hi,Pi,h,'cubic');
    for(i=1:1:N)
        y(i) = interp1(hi,yi(:,i),h,'cubic');
    end
    V = interp1(hi,Vi,h,'cubic');
    rf = interp1(hi,rfi,h,'cubic');
    rg = interp1(hi,rgi,h,'cubic');
end

htoler = 1;     % J/kg
Pinc = 10;      % Use a small increment for P derivatives, Pa.
Plast = P;      % Save last value for convergence on P.

imax = 10;
for(i=1:1:imax)
    P
    % See how we are doing on P.
    [q V y x rg rf] = Flash_zTP(z,T,P,V,y,rf,rg,1);
    hf = h_crT(x,rf,T);
    hg = h_crT(y,rg,T);
    hest = (1-q)*hf + q*hg;
    f = hest - h;

    if((abs(f) < htoler)&&((abs(P-Plast)/P) < toler))
        break
    end
    Plast = P;
    
    % Use NR to adjust pressure to get h correct.
    % Use a central difference for the numerical derivative.
    Plow = P - Pinc;
    [q V y x rg rf] = Flash_zTP(z,T,Plow,V,y,rf,rg,1);
    hf = h_crT(x,rf,T);
    hg = h_crT(y,rg,T);
    hlow = (1-q)*hf + q*hg;
    
    Phigh = P + Pinc;
    [q V y x rg rf] = Flash_zTP(z,T,Phigh,V,y,rf,rg,1);
    hf = h_crT(x,rf,T);
    hg = h_crT(y,rg,T);
    hhigh = (1-q)*hf + q*hg;
    
    dfdP = ((hhigh-h)-(hlow-h))/(2*Pinc);
    P = P - f/dfdP;
end
iP = i
if(i == imax)
    disp('Fell off the end of NR loop for P in Flash_zhT')
end
P = P;
q = q;
V = V;
y = y;
x = x;
rg = rg;
rf = rf;