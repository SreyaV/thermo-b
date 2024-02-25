function [P rg rf x] = Dew_cT(y,T,varargin)
% Return the dew-point pressure P (Pa), molefractions x, and densities (kg/m3) 
% for given vapor molefractions y and temperature T (K).
% C.F. Edwards, 2-13-10

global toler

% Use information supplied if available.
switch nargin
    case 2
        [Pstart rvstart rlstart xstart] = Dew_Start_cT(y,T);
    case 6
        Pstart  = varargin{1};
        xstart  = varargin{2};
        rlstart = varargin{3};
        rvstart = varargin{4};
    otherwise
        disp('Incorrect number of arguments in Dew_cT')
        return
end

Pinc = 10;      % Use a small increment for P derivatives.
P = Pstart;     % Initialize the pressure.
Plast = P;      % Save last value for convergence on P.
x = xstart;
rl = rlstart;
rv = rvstart;

% Use the starting composition to determine which species to leave for use
% in closure.  This should be the one with the highest mole fraction since
% it will leave the most head room against NR trying to exceed unity for
% any species.
[maxsp isp] = max(xstart);

imax = 100;
for i=1:1:imax
    % Find the target chemical potential for isp in the vapor.
    rv = rv_cTP(y,T,P,rv);
    muvisp = mui_icrT(isp,y,rv,T);

    % Find the matching liquid compositions for all but isp.
    x = x_iyTPr(isp,y,T,P,rl,rv,x);
    
    % See how we are doing on isp.
    rl = rl_cTP(x,T,P,rl);
    mulisp = mui_icrT(isp,x,rl,T);
    fisp = mulisp - muvisp;

    if((abs(fisp/muvisp) < toler)&&((abs(P-Plast)/P) < toler))
        rg = rv;
        rf = rl;
        x = real(x);
        return
    end
    Plast = P;
    
    % Use NR to adjust pressure to get isp correct.
    % Choose the order of the derivative by setting the next line.
    order = 2;
    if(order == 2)
        % Use a second-order central difference for the numerical derivative.
        % Go low side.
        rv = rv_cTP(y,T,P-Pinc,rv);
        muvlow = mui_icrT(isp,y,rv,T);
        x = x_iyTPr(isp,y,T,P-Pinc,rl,rv,x);
        rl = rl_cTP(x,T,P-Pinc,rl);
        mullow = mui_icrT(isp,x,rl,T);
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvhigh = mui_icrT(isp,y,rv,T);
        x = x_iyTPr(isp,y,T,P+Pinc,rl,rv,x);
        rl = rl_cTP(x,T,P+Pinc,rl);
        mulhigh = mui_icrT(isp,x,rl,T);
        % Form the derivative.
        dfispdP = ((mulhigh-muvhigh)-(mullow-muvlow))/(2*Pinc);
    else
        % Use a first-order forward difference for the numerical derivative.
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvhigh = mui_icrT(isp,y,rv,T);
        x = x_iyTPr(isp,y,T,P+Pinc,rl,rv,x);
        rl = rl_cTP(x,T,P+Pinc,rl);
        mulhigh = mui_icrT(isp,x,rl,T);
        % Form the derivative.
        dfispdP = ((mulhigh-muvhigh)-(mulisp-muvisp))/(Pinc);
    end

    % Get the new pressure estimate.
    P = P - fisp/dfispdP;
end

disp('Fell off the end of NR loop for P in Dew_Point_cT')
