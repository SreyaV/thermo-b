function [T rg rf x] = Dew_cP(y,P,varargin)
% Return the dew-point pressure P (Pa), molefractions x, and densities (kg/m3) 
% for given vapor molefractions y and temperature T (K).
% C.F. Edwards, 2-16-10

global toler

% Use information supplied if available.
switch nargin
    case 2
        % The fact that most saturation curves are nearly linear on a plot of
        % lnP vs. 1/T can be used to estimate the initial temperature.
        % Use the inflection point as the upper end.
        [Tinfl rinfl] = Pr_Inflection_c(y);
        Pinfl = P_crT(y,rinfl,Tinfl);
        % This first T is low enough to use an ideal solution.
        Tlow = 70;
        Plow = Ideal_Dew_cT(y,Tlow);
        % Change to the following line if that does not work.
        % Plow = Dew_cT(y,Tlow);
        Tstart = 1/(1/Tinfl - (1/Tlow - 1/Tinfl)*log(P/Pinfl)/log(Pinfl/Plow));
        [Pstart rvstart rlstart xstart] = Dew_cT(y,Tstart);
    case 6
        Tstart  = varargin{1};
        xstart  = varargin{2};
        rlstart = varargin{3};
        rvstart = varargin{4};
    otherwise
        disp('Incorrect number of arguments in Dew_Point_cP')
        return
end

Tinc = Tstart*0.001;    % Use a small increment for derivatives.
T = Tstart;             % Initialize.
Tlast = T;              % Save last value for convergence.
x = xstart;
rl = rlstart;
rv = rvstart;
    
imax = 100;
for i=1:1:imax
    T
    % Find the pressure for this temperature estimate.
    [Pest rv rl x] = Dew_cT(y,T,P,x,rl,rv);
    
    % See how we are doing.
    f = P - Pest;

    if((abs(f/P) < toler)&&((abs(T-Tlast)/T) < toler))
        break
    end
    Tlast = T;
    
    % Use NR to adjust temperature to get P correct.
    % Use a central difference for the numerical derivative.
    [Phigh rv rl x] = Dew_cT(y,T+Tinc,Pest,x,rl,rv);
    
    dfdT = ((P-Phigh)-(P-Pest))/Tinc;
    T = T - f/dfdT;
end

if(i == imax)
    disp('Fell off the end of NR loop for T in Dew_cP')
end
rg = rv;
rf = rl;
