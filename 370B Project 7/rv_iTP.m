function r = rv_iTP(ispecies,T,P,varargin)
% Return the vapor density (kg/m3) of species i for any given T (K) and P (Pa).
% C.F. Edwards, 2-3-10

global toler
global Tcrit_i rcrit_i rgtrip_i Pcrit_i

% Shorten the names.
Tcrit = Tcrit_i(ispecies);
rcrit = rcrit_i(ispecies);
Pcrit = Pcrit_i(ispecies);
rgtrip = rgtrip_i(ispecies);

% Set the tolerance for the convergence of the density.  Use the critical
% density as the representative value.
rtoler = sqrt(eps(rcrit));
ftoler = 1; % Pa

% Check for near critical.
if((abs(T-Tcrit)/T < toler)&&(abs(P-Pcrit)/P < toler))
    r = rcrit;
    return
end

% Set a maximum step length.  The intent here is to avoid excessively large
% steps that are clearly non-physical.  Try a fraction of the half-distance
% across the dome.
drmax = (rcrit-rgtrip)/2;

% Use a supplied starting point if available.  If not, use the triple line
% density on the vapor side (or some fraction thereof).
if(nargin == 3)    
    rstart = rgtrip/10;
else
    rstart = varargin{1};
end

% Use Newton-Raphson to find the density.  We will approach this as a
% zero-crossing problem for P using the approach outlined in Chapter 2 of
% Numerical Methods for Unconstrained Optimization and Nonlinear Equations
% by J.E. Dennis and R.B. Schnabel, Prentice-Hall, 1983.

% Set the parameter for required amount of decrease in the function.
alpha = 1e-4;

% Set how many times we will try the N-R iteration before complaining.
imax = 20;

% Set the initial values of the density and find the function and its
% derivative there.
r = rstart;
f = P_irT(ispecies,r,T)-P;

% For the first pass, set the "last" value to the current value.  This will
% cause an immediate return only if the function is already within toler of
% a zero.
rlast = r;

% Start the N-R iteration loop.
for i=1:1:imax
    % Test to see if this is the answer.
    f_small_enough = (abs(f) < ftoler);
    r_close_enough    = (abs(r-rlast)/r < rtoler);
    if(f_small_enough && r_close_enough)
%         Iterations = i
        r = r;
        return
    end
    % Not close enough.  Try a full Newton step (drN).
    % This is based on a linear projection using the derivative.
    dfdr = dPdr_irT(ispecies,r,T);
    % Find the Newton displacement.  Note that the step depends on the
    % ratio of the derivative to the magnitude of the second derivative.
    % This way you always move in the direction that reduces the value
    % of the function (not its derivative).
    drN = -f/abs(dfdr);
    % Limit the step size if needed.  Keep the direction.
    if(abs(drN) > drmax)
%         disp('Hit max step limiter in rv_iTP')
        drN = drmax*drN/abs(drN);
    end
    % Test to see if this step is acceptable.  If so, continue around the
    % N-R loop.  If not, use a backtracking method to give an accetable
    % step in the Newton direction.
    [r f] = Linesearch(ispecies,T,P,r,f,dfdr,drN,alpha);
    rlast = r;
end

% Were we successful?
if(i == imax)
    % Nope.  Sorry!
    disp('Zero not found in rv_iTP')
    r = 0;
    return
else
    % Yep!  Here it is.
    r = r;
    return
end


function [r f] = Linesearch(ispecies,T,P,r,f,dfdr,drN,alpha)
% The criteria for an acceptable step is that the magnitude of the function 
% decrease by some fraction alpha of what you would get from a full Newton
% step (as evaluated by the derivative using the current position).

% See if the Newton step satisfies the decrease criterion.
rN = r + drN;
fN = P_irT(ispecies,rN,T)-P;
df_big_enough = ((abs(f)-abs(fN)) > -alpha*dfdr*drN);
if(df_big_enough)
    % Accept the full Newton step.
    r = rN;
    f = fN;
    return
end

% Do not accept the full Newton step.
% disp('Step not accepted')

% Backtrack to a smaller but acceptable step.

% Set how many times to try before complaining.
jmax = 20;

% Mark the previous attempt as being the Newton step.
lambda = 1;
drL = drN;

for j=1:1:jmax
    % Try a half step.
    lambda = lambda/2;
    drL = lambda*drL;                   % Lambda step
    rL  = r + drL;                      % Lambda location
    fL  = P_irT(ispecies,rL,T) - P;     % Function at Lambda
    
    % See if the lambda step satisfies the decrease criterion.
    df_big_enough = ((abs(f)-abs(fL)) > -alpha*dfdr*drL);
    if(df_big_enough)
        % Accept the lambda step.
%         Iterations = j
        r = rL;
        f = fL;
        return
    end
end
if(j == jmax)
%     disp('Backtrack not successfull in rv_iTP')
end
