function [T rf rg y] = Bubble_cP(x,P,varargin)
% Return the bubble-point temperature T (K), liquid and vapor densities (kg/m3),
% and vapor molefractions y for given liquid molefractions x and pressure P (Pa).
% C.F. Edwards, 2-16-10

global toler

% Use information supplied if available.
switch nargin
    case 2
        % The fact that most saturation curves are nearly linear on a plot of
        % lnP vs. 1/T can be used to estimate the initial temperature.
        [Tinfl rinfl] = Pr_Inflection_c(x);
        Pinfl = P_crT(x,rinfl,Tinfl);
        Tlow = 70;
        % This first T is low enough to use an ideal solution.
        Plow = Ideal_Bubble_cT(x,Tlow);
        % Change to the following line if that does not work.
        % Plow = Bubble_Point_cT(x,Tlow);
        Tstart = 1/(1/Tinfl - (1/Tlow - 1/Tinfl)*log(P/Pinfl)/log(Pinfl/Plow));
        [Pstart rlstart rvstart ystart] = Bubble_cT(x,Tstart);
    case 6
        Tstart  = varargin{1};
        ystart  = varargin{2};
        rlstart = varargin{3};
        rvstart = varargin{4};
    otherwise
        disp('Incorrect number of arguments in Bubble_cP')
        return
end

Tinc = Tstart*0.001;    % Use a small increment for derivatives.
T = Tstart;             % Initialize.
Tlast = T;              % Save last value for convergence.
y = ystart;
rl = rlstart;
rv = rvstart;

imax = 100;
for i=1:1:imax
    % Find the pressure for this temperature estimate.
    [Pest rl rv y] = Bubble_cT(x,T,P,y,rl,rv);
    
    % See how we are doing.
    f = P - Pest;

    if((abs(f/P) < toler)&&((abs(T-Tlast)/T) < toler))
        break
    end
    Tlast = T;
    
    % Use NR to adjust temperature to get P correct.
    % Use a central difference for the numerical derivative.
    [Phigh rl rv y] = Bubble_cT(x,T+Tinc,Pest,y,rl,rv);
    
    dfdT = ((P-Phigh)-(P-Pest))/Tinc;
    T = T - f/dfdT;
end

if(i == imax)
    disp('Fell off the end of NR loop for T in Bubble_cP')
end
rf = rl;
rg = rv;
