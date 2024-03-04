function [P rf rg y] = Ideal_Bubble_cT(x,T,varargin)
% Return the bubble-point pressure P (Pa), molefractions y, and densities (kg/m3) 
% for given liquid molefractions x and temperature T (K) as calculated by
% assuming that the vapor complement can be treated as an ideal solution.
% C.F. Edwards, 2-13-10

% No starting values are required.  The purpose of this function is to
% generate them.

global toler
global N2 O2 Ar
global Ru M_i

% Use information supplied if available.
switch nargin
    case 2
        % Get the inflection-point properties for the liquid composition.
        % Doing this once will save time on rl_cTP calls.
        [Tinfl_x rinfl_x] = Pr_Inflection_c(x);
    case 4
        Tinfl_x = varargin{1};
        rinfl_x = varargin{2};
    otherwise
        disp('Incorrect number of arguments in Ideal_Dew_cT')
        return
end
Pinfl_x = P_crT(x,rinfl_x,Tinfl_x);

% Find out how many components.
N = length(x);  
if ~((N == 2)||(N == 3))
    disp('Number of components must be 2 or 3 in Ideal_Bubble_cT')
    return
end

% Preallocate some storage.
y     = zeros(1,N);
ylow  = zeros(1,N);
yhigh = zeros(1,N);
Pgsi  = zeros(1,N);

% Use a small increment for P derivatives.
Pinc = 10;  % Pascals

% Find the highest pressure at which you can make an ideal solution.
for i=1:1:length(x)
    rgs = Vapor_Spinodal_iT(i,T);
    if(rgs == 0)
        Pgsi(i) = Pinfl_x;
    else
        Pgsi(i) = P_irT(i,rgs,T);
    end
end
% The smallest of these values sets the upper bound.
Pmax = 0.99*min(Pgsi);

% Find the lowest pressure at which you can evaluate the liquid mixture
% chemical potentials.
rfs = Liquid_Spinodal_cT(x,T);
Pfs = P_crT(x,rfs,T);
% Watch out for bad spinodals!
if(Pfs < 0)
    Pmin = 1;
else
    Pmin = Pfs;
end

% Check to see if an ideal solution is possible.  If not, return zeros.
if(Pmax <= Pmin)
    disp('Pressures do not overlap in Ideal_Bubble_cT')
    P = 0; rf = 0; rg = 0; y = zeros(1,N);
    return
end

% Choose a starting value in this range.
Pstart = 0.5*(Pmax-Pmin) + Pmin;

% Set at starting pressure and liquid density.
P = Pstart;
rl = rl_cTP(x,T,P,Tinfl_x,rinfl_x);
rvN2 = rv_iTP(N2,T,P);
rvO2 = rv_iTP(O2,T,P);
if(N == 3)
    % Argon is included.
    rvAr = rv_iTP(Ar,T,P);
end
    
imax = 100;
Plast = P;
for i=1:1:imax 
    % Find the target chemical potentials in the liquid.
    % Note that these are the actual, mixture-model values, not an ideal
    % approximation.
    rl = rl_cTP(x,T,P,Tinfl_x,rinfl_x,rl);
    mulO2 = mui_icrT(O2,x,rl,T);
    mulN2 = mui_icrT(N2,x,rl,T);
    
    % Find the chemical potentials of the pure vapor components at this
    % thermal state (T,P).
    rvN2 = rv_iTP(N2,T,P,rvN2);
    muvN2neat = mu_irT(N2,rvN2,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(N2,[1 0 0],rvN2,T), but this is faster.
    rvO2 = rv_iTP(O2,T,P,rvO2);
    muvO2neat = mu_irT(O2,rvO2,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(O2,[0 1 0],rvO2,T), but this is faster.
    
    % Set the mole fractions as per an ideal solution.
    % Note that this is an ideal solution only for the complementary phase,
    % since the known phase is, known, we don't need to approximate it.
    y(O2) = exp((mulO2 - muvO2neat)/Ru/T);
    y(N2) = exp((mulN2 - muvN2neat)/Ru/T);
    
    if(N == 3)
        % Argon is included.
        mulAr = mui_icrT(Ar,x,rl,T);
        rvAr = rv_iTP(Ar,T,P,rvAr);
        muvArneat = mu_irT(Ar,rvAr,T);
        % Note that you could also get the value above by using the mixture
        % function mui_icrT(Ar,[0 0 1],rvAr,T), but this is faster.
        y(Ar) = exp((mulAr - muvArneat)/Ru/T);
    end
    
    % Make the composite vapor density via Amagat.
    mN2 = y(N2)*M_i(N2);
    vN2 = mN2/rvN2;
    mO2 = y(O2)*M_i(O2);
    vO2 = mO2/rvO2;
    vmix = vN2 + vO2;
    mmix = mN2 + mO2;
    if(N == 3)
        mAr = y(Ar)*M_i(Ar);
        vAr = mAr/rvAr;
        vmix = vmix + vAr;
        mmix = mmix + mAr;
    end
    rvmix = mmix/vmix;
        
    % Check on the sum to see if it adds to unity.
    ysum = sum(y);
    f = ysum - 1;
    
    % Are we done?
    if((abs(f) < 1e-3)&&(((abs(P-Plast)/P) < toler)||(abs(P-Plast) < 1)))
        rg = rvmix;
        rf = rl;
        y = real(y/ysum);
        return
    end
    Plast = P;
    
    % Use NR to adjust pressure to get the mole fraction sum correct.
    % Choose the order of the derivative by setting the next line.
    order = 2;
    if(order == 2)
        % Use a second-order central difference for the numerical derivative.
        % Go low side.
        rl = rl_cTP(x,T,P-Pinc,Tinfl_x,rinfl_x,rl);
        mulO2 = mui_icrT(O2,x,rl,T);
        mulN2 = mui_icrT(N2,x,rl,T);
        rvN2 = rv_iTP(N2,T,P-Pinc,rvN2);
        muvN2neat = mu_irT(N2,rvN2,T);
        rvO2 = rv_iTP(O2,T,P-Pinc,rvO2);
        muvO2neat = mu_irT(O2,rvO2,T);
        ylow(O2) = exp((mulO2 - muvO2neat)/Ru/T);
        ylow(N2) = exp((mulN2 - muvN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            mulAr = mui_icrT(Ar,x,rl,T);
            rvAr = rv_iTP(Ar,T,P-Pinc,rvAr);
            muvArneat = mu_irT(Ar,rvAr,T);
            ylow(Ar) = exp((mulAr - muvArneat)/Ru/T);
        end
        ysumlow = sum(ylow);
        flow = ysumlow - 1;
        % Go high side.
        rl = rl_cTP(x,T,P+Pinc,Tinfl_x,rinfl_x,rl);
        mulO2 = mui_icrT(O2,x,rl,T);
        mulN2 = mui_icrT(N2,x,rl,T);
        rvN2 = rv_iTP(N2,T,P+Pinc,rvN2);
        muvN2neat = mu_irT(N2,rvN2,T);
        rvO2 = rv_iTP(O2,T,P+Pinc,rvO2);
        muvO2neat = mu_irT(O2,rvO2,T);
        yhigh(O2) = exp((mulO2 - muvO2neat)/Ru/T);
        yhigh(N2) = exp((mulN2 - muvN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            mulAr = mui_icrT(Ar,x,rl,T);
            rvAr = rv_iTP(Ar,T,P+Pinc,rvAr);
            muvArneat = mu_irT(Ar,rvAr,T);
            yhigh(Ar) = exp((mulAr - muvArneat)/Ru/T);
        end
        ysumhigh = sum(yhigh);
        fhigh = ysumhigh - 1;
        % Form the derivative.
        dfdP = (fhigh - flow)/(2*Pinc);
    else
        % Use a first-order forward difference for the numerical derivative.
        % Go high side.
        rl = rl_cTP(x,T,P+Pinc,Tinfl_x,rinfl_x,rl);
        mulO2 = mui_icrT(O2,x,rl,T);
        mulN2 = mui_icrT(N2,x,rl,T);
        rvN2 = rv_iTP(N2,T,P+Pinc,rvN2);
        muvN2neat = mu_irT(N2,rvN2,T);
        rvO2 = rv_iTP(O2,T,P+Pinc,rvO2);
        muvO2neat = mu_irT(O2,rvO2,T);
        yhigh(O2) = exp((mulO2 - muvO2neat)/Ru/T);
        yhigh(N2) = exp((mulN2 - muvN2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            mulAr = mui_icrT(Ar,x,rl,T);
            rvAr = rv_iTP(Ar,T,P+Pinc,rvAr);
            muvArneat = mu_irT(Ar,rvAr,T);
            yhigh(Ar) = exp((mulAr - muvArneat)/Ru/T);
        end
        ysumhigh = sum(yhigh);
        fhigh = ysumhigh - 1;
        % Form the derivative.
        dfdP = ((fhigh)-(f))/(Pinc);
    end

    % Get the new pressure estimate.
    dP = -f/dfdP;
    % Limit the step size.
    limit = (Pmax+Pmin)/4;
    if(abs(dP) > limit)
%         disp('Hit step limiter in Ideal_Bubble_cT') 
        dP = sign(dP)*limit;
    end
    P = P + dP;
    % Watch for out of bounds.
    if(P < Pmin)
        P = Pmin;
    end
    if(P > Pmax)
        P = Pmax;
    end
end

P = 0; rf = 0; rg = 0; y = zeros(1,N);
disp('Fell off the end of NR loop for P in Ideal_Bubble_cT')
