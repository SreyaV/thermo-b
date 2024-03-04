function [P rg rf x] = Ideal_Dew_cT(y,T,varargin)
% Return the dew-point pressure P (Pa), molefractions x, and densities (kg/m3) 
% for given vapor molefractions y and temperature T (K) as calculated by
% assuming that the liquid complement can be treated as an ideal solution.
% C.F. Edwards, 2-13-10

% No starting values are required.  The purpose of this function is to
% generate them.

global toler
global N2 O2 Ar
global Ru M_i

% Use information supplied if available.
switch nargin
    case 2
        % Get the inflection-point properties for the vapor composition.  Doing this
        % once will save time on rv_cTP calls.
        [Tinfl_y rinfl_y] = Pr_Inflection_c(y);
    case 4
        Tinfl_y = varargin{1};
        rinfl_y = varargin{2};
    otherwise
        disp('Incorrect number of arguments in Ideal_Dew_cT')
        return
end
Pinfl_y = P_crT(y,rinfl_y,Tinfl_y);

% Find out how many components.
N = length(y);  
if ~((N == 2)||(N == 3))
    disp('Number of components must be 2 or 3 in Ideal_Dew_cT')
    return
end

% Preallocate some storage.
x     = zeros(1,N);
xlow  = zeros(1,N);
xhigh = zeros(1,N);
Pfsi  = zeros(1,N);

% Use a small increment for P derivatives.
Pinc = 10;  % Pascals

% Find the lowest pressure at which you can make an ideal solution.
for i=1:1:length(y)
    rfs = Liquid_Spinodal_iT(i,T);
    Pfsi(i) = P_irT(i,rfs,T);
end
% The largest of these sets the lower bound.  If all are negative, use
% a small value (like 1 Pa).
Pmin = 1.01*max([Pfsi 1]);

% Find the highest pressure at which you can evaluate the vapor mixture
% chemical potentials.
rgs = Vapor_Spinodal_cT(y,T);
Pgs = P_crT(y,rgs,T);
% Watch out for bad spinodals!
if(Pgs > Pinfl_y)
    Pmax = Pinfl_y;
else
    Pmax = Pgs;
end

% Check to see if an ideal solution is possible.  If not, return zeros.
if(Pmax < Pmin)
    disp('Pressures do not overlap in Ideal_Dew_cT')
    P = 0; rg = 0; rf = 0; x = [0 0 0];
    return
end

% Choose a starting value in this range.
Pstart = 0.5*(Pmax-Pmin) + Pmin;

% Set at starting pressure and vapor density.
P = Pstart;
rv = rv_cTP(y,T,P,Tinfl_y,rinfl_y);
rlN2 = rl_iTP(N2,T,P);
rlO2 = rl_iTP(O2,T,P);
if(N == 3)
    % Argon is included.
    rlAr = rl_iTP(Ar,T,P);
end
    
imax = 50;
Plast = P;
for i=1:1:imax
    % Find the target chemical potentials in the vapor.
    % Note that these are the actual, mixture-model values, not an ideal
    % approximation.
    rv = rv_cTP(y,T,P,Tinfl_y,rinfl_y,rv);
    muvO2 = mui_icrT(O2,y,rv,T);
    muvN2 = mui_icrT(N2,y,rv,T);
    
    % Find the chemical potentials of the pure liquid components.
    rlN2 = rl_iTP(N2,T,P,rlN2);
    mulN2neat = mu_irT(N2,rlN2,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(N2,[1 0 0],rvN2,T), but this is faster.
    rlO2 = rl_iTP(O2,T,P,rlO2);
    mulO2neat = mu_irT(O2,rlO2,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(O2,[0 1 0],rvO2,T), but this is faster.
    
    % Set the mole fractions as per an ideal solution.
    x(O2) = exp((muvO2 - mulO2neat)/Ru/T);
    x(N2) = exp((muvN2 - mulN2neat)/Ru/T);
    
    if(N == 3)
        % Argon is included.
        muvAr = mui_icrT(Ar,y,rv,T);
        rlAr = rl_iTP(Ar,T,P,rlAr);
        mulArneat = mu_irT(Ar,rlAr,T);
        % Note that you could also get the value above by using the mixture
        % function mui_icrT(Ar,[0 0 1],rvAr,T), but this is faster.
        x(Ar) = exp((muvAr - mulArneat)/Ru/T);
    end

    % Make the composite liquid density via Amagat.
    mN2 = x(N2)*M_i(N2);
    vN2 = mN2/rlN2;
    mO2 = x(O2)*M_i(O2);
    vO2 = mO2/rlO2;
    vmix = vN2 + vO2;
    mmix = mN2 + mO2;
    if(N == 3)
        mAr = x(Ar)*M_i(Ar);
        vAr = mAr/rlAr;
        vmix = vmix + vAr;
        mmix = mmix + mAr;
    end
    rlmix = mmix/vmix;
    
    % Check on the sum to see if it adds to unity.
    xsum = sum(x);
    f = xsum - 1;
    
    % Are we done?
    if((abs(f) < toler)&&((abs(P-Plast)/P) < toler))
        rg = rv;
        rf = rlmix;
        x = real(x/xsum);
        return
    end
    Plast = P;
    
    % Use NR to adjust pressure to get the mole fraction sum correct.
    % Choose the order of the derivative by setting the next line.
    order = 2;
    if(order == 2)
        % Use a second-order central difference for the numerical derivative.
        % Go low side.
        rv = rv_cTP(y,T,P-Pinc,Tinfl_y,rinfl_y,rv);
        muvO2 = mui_icrT(O2,y,rv,T);
        muvN2 = mui_icrT(N2,y,rv,T);
        rl = rl_iTP(N2,T,P-Pinc,rlN2);
        mulN2neat = mu_irT(N2,rl,T);
        rl = rl_iTP(O2,T,P-Pinc,rlO2);
        mulO2neat = mu_irT(O2,rl,T);
        xlow(O2) = exp((muvO2 - mulO2neat)/Ru/T);
        xlow(N2) = exp((muvN2 - mulN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P-Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xlow(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumlow = sum(xlow);
        flow = xsumlow - 1;
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,Tinfl_y,rinfl_y,rv);
        muvO2 = mui_icrT(O2,y,rv,T);
        muvN2 = mui_icrT(N2,y,rv,T);
        rl = rl_iTP(N2,T,P+Pinc,rlN2);
        mulN2neat = mu_irT(N2,rl,T);
        rl = rl_iTP(O2,T,P+Pinc,rlO2);
        mulO2neat = mu_irT(O2,rl,T);
        xhigh(O2) = exp((muvO2 - mulO2neat)/Ru/T);
        xhigh(N2) = exp((muvN2 - mulN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P+Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xhigh(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumhigh = sum(xhigh);
        fhigh = xsumhigh - 1;
        % Form the derivative.
        dfdP = ((fhigh)-(flow))/(2*Pinc);
    else
        % Use a first-order forward difference for the numerical derivative.
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,Tinfl_y,rinfl_y,rv);
        muvO2 = mui_icrT(O2,y,rv,T);
        muvN2 = mui_icrT(N2,y,rv,T);
        rl = rl_iTP(N2,T,P+Pinc,rlN2);
        mulN2neat = mu_irT(N2,rl,T);
        rl = rl_iTP(O2,T,P+Pinc,rlO2);
        mulO2neat = mu_irT(O2,rl,T);
        xhigh(O2) = exp((muvO2 - mulO2neat)/Ru/T);
        xhigh(N2) = exp((muvN2 - mulN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P+Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xhigh(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumhigh = sum(xhigh);
        fhigh = xsumhigh - 1;
        % Form the derivative.
        dfdP = ((fhigh)-(f))/(Pinc);
    end

    % Get the new pressure estimate.
    dP = -f/dfdP;
    % Limit the step size.
%     limit = 100e3;
%     if(abs(dP) > limit)
%         disp('Hit step limiter in Ideal_Dew_cT') 
%         dP = sign(dP)*limit;
%     end
    P = P + dP;
    % Watch for out of bounds.
    if(P < Pmin)
        P = Pmin;
    end
    if(P > Pmax)
        P = Pmax;
    end
end

P = 0; rg = 0;  rf = 0; x = zeros(1,N);
disp('Fell off the end of NR loop for P in Ideal_Dew_cT')
