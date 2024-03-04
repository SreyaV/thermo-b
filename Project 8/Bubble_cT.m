function [P rf rg y] = Bubble_cT(x,T,varargin)
% Return the bubble-point pressure P (Pa), liquid and vapor densities (kg/m3),
% and vapor molefractions y for given liquid molefractions x and temperature T (K).
% C.F. Edwards, 2-13-10

global toler

% Use information supplied if available.
switch nargin
    case 2
        [Pstart rlstart rvstart ystart] = Bubble_Start_cT(x,T);
    case 6
        Pstart  = varargin{1};
        ystart  = varargin{2};
        rlstart = varargin{3};
        rvstart = varargin{4};
    otherwise
        disp('Incorrect number of arguments in Bubble_cT')
        return
end

Pinc = 10;      % Use a small increment for P derivatives.
P = Pstart;     % Initialize the pressure.
Plast = P;      % Save last value for convergence on P.
y = ystart;
rl = rlstart;
rv = rvstart;

% Use the starting composition to determine which species to leave for use
% in closure.  This should be the one with the highest mole fraction since
% it will leave the most head room against NR trying to exceed unity for
% any species.
[maxsp isp] = max(ystart);

imax = 100;
for i=1:1:imax
    % Find the target chemical potential for isp in the liquid.
    rl = rl_cTP(x,T,P,rl);
    mulisp = mui_icrT(isp,x,rl,T);

    % Find the matching vapor compositions for all but isp.
    y = y_ixTPr(isp,x,T,P,rl,rv,y);
    
    % See how we are doing on isp.
    rv = rv_cTP(y,T,P,rv);
    muvisp = mui_icrT(isp,y,rv,T);
    fisp = muvisp - mulisp;

    % It is hard to set the convergence criterion for mu.
    % Find the mu error that corresponds to a T diff for the vapor.
    dT = toler*300;
    dmuv = abs(muvisp-mui_icrT(isp,y,rv,T+dT));
    if((abs(fisp) < dmuv)&&((abs(P-Plast)/P) < toler))
        rf = rl;
        rg = rv;
        y = real(y);
        return
    end
    Plast = P;
    
    % Use NR to adjust pressure to get isp correct.
    % Choose the order of the derivative by setting the next line.
    order = 2;
    if(order == 2)
        % Use a second-order central difference for the numerical derivative.
        % Go low side.
        rl = rl_cTP(x,T,P-Pinc,rl);
        mullow = mui_icrT(isp,x,rl,T);
        y = y_ixTPr(isp,x,T,P-Pinc,rl,rv,y);
        rv = rv_cTP(y,T,P-Pinc,rv);
        muvlow = mui_icrT(isp,y,rv,T);
        % Go high side.
        rl = rl_cTP(x,T,P+Pinc,rl);
        mulhigh = mui_icrT(isp,x,rl,T);
        y = y_ixTPr(isp,x,T,P+Pinc,rl,rv,y);
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvhigh = mui_icrT(isp,y,rv,T);
        % Form the derivative.
        dfispdP = ((muvhigh-mulhigh)-(muvlow-mullow))/(2*Pinc);
    else
        % Use a first-order forward difference for the numerical derivative.
        % Go high side.
        rl = rl_cTP(x,T,P+Pinc,rl);
        mulhigh = mui_icrT(isp,x,rl,T);
        y = y_ixTPr(isp,x,T,P+Pinc,rl,rv,y);
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvhigh = mui_icrT(isp,y,rv,T);
        % Form the derivative.
        dfispdP = ((muvhigh-mulhigh)-(muvisp-mulisp))/(Pinc);
    end
    
    % Get the new pressure estimate.
    P = P - fisp/dfispdP;
end

disp('Fell off the end of NR loop for P in Bubble_Point_cT')

