function rfs = Liquid_Spinodal_cT(c,T,varargin)
% Return the density of the liquid-side spinodal line rfs (kg/m3) for
% composition c at temperature T (K).

% If there are additional arguments specified, they are:
% a starting density (kg/m3), or
% Tinfl (K) and rinfl (kg/m3), or
% Tinfl (K), rinfl (kg/m3), and a starting density (kg/m3).

% If the temperature is above inflection or the spinodal cannot be found,
% return a value of zero for the density and complain.
% C.F. Edwards, 2-16-10

global toler

% Use supplied information if available.
switch nargin
    case 5
        Tinfl  = varargin{1};
        rinfl  = varargin{2};
        rstart = varargin{3};
    case 4
        Tinfl  = varargin{1};
        rinfl  = varargin{2};
        rstart = 2*rinfl;
    case 3
        rstart = varargin{1};
        [Tinfl rinfl] = Pr_Inflection_c(c);
    otherwise
        [Tinfl rinfl] = Pr_Inflection_c(c);
        rstart = 2*rinfl;
end

% Set the tolerance for the convergence of the density.  Use the inflection
% density as the representative value.
rtoler = sqrt(eps(rinfl));
ftoler = 1; % Pa/kg/m3

% Check for near inflection.
if(abs(T-Tinfl)/T < toler)
    rfs = rinfl;
    return
end

% Check for T > Tinfl.
if(T > Tinfl)
    disp('Temperature above inflection in Liquid_Spinodal_cT')
    rfs = 0;
    return
end

% Set a maximum step length.  The intent here is to avoid excessively large
% steps that are clearly non-physical.  Try a fraction of the half-distance
% across the dome.
drmax = (rinfl)/4;

% Use Newton-Raphson to find the spinodal.  We will approach this as a
% minimization problem for P using the approach outlined in Chapter 2 of
% Numerical Methods for Unconstrained Optimization and Nonlinear Equations
% by J.E. Dennis and R.B. Schnabel, Prentice-Hall, 1983.

% Set the parameter for required amount of decrease in the function.
alpha = 1e-4;

% Set how many times we will try the N-R iteration before complaining.
imax = 20;

% Set the initial values of the density and find the function and its
% derivative there.
r = rstart;
dfdr = dPdr_crT(c,r,T);

% For the first pass, set the "last" value to the current value.  This will
% cause an immediate return only if the function is already within toler of
% a zero.
rlast = r;

% Start the N-R iteration loop.
for i=1:1:imax
    % Test to see if this is the answer.
    dfdr_small_enough = (abs(dfdr) < ftoler);
    r_close_enough    = (abs(r-rlast)/r < rtoler);
    if(dfdr_small_enough && r_close_enough)
%         Iterations = i
        rfs = r;
        return
    end
    % Not close enough.  Try a full Newton step (drN).
    % This is based on a linear projection using the derivative.
    d2fdr2 = d2Pdr2_crT(c,r,T);
    % Find the Newton displacement.  Note that the step depends on the
    % ratio of the derivative to the magnitude of the second derivative.
    % This way you always move in the direction that reduces the value
    % of the function (not its derivative).
    drN = -dfdr/abs(d2fdr2);
    % Limit the step size if needed.  Keep the direction.
    if(abs(drN) > drmax)
%         disp('Hit max step limiter in Liquid_Spinodal_iT')
        drN = drmax*drN/abs(drN);
    end
    % Test to see if this step is acceptable.  If so, continue around the
    % N-R loop.  If not, use a backtracking method to give an accetable
    % step in the Newton direction.
    [r dfdr] = Linesearch(c,T,r,dfdr,drN,alpha);
    rlast = r;
end

% Were we successful?
if(i == imax)
    % Nope.  Sorry!
    disp('Spinodal not found in Liquid_Spinodal_cT')
    rfs = 0;
    return
else
    % Yep!  Here it is.
    rfs = r;
    return
end


function [r dfdr] = Linesearch(c,T,r,dfdr,drN,alpha)
% The criteria for an acceptable step is that the magnitude of the function 
% decrease by some fraction alpha of what you would get from a full Newton
% step (as evaluated by the derivative using the current position).

% See if the Newton step satisfies the decrease criterion.
f = P_crT(c,r,T);
rN = r + drN;
fN = P_crT(c,rN,T);
df_big_enough = ((f-fN) > -alpha*dfdr*drN);
if(df_big_enough)
    % Accept the full Newton step.
    r = rN;
    dfdr = dPdr_crT(c,rN,T);
    return
end

% Do not accept the full Newton step.
% disp('Step not accepted')
% Backtrack to a smaller but acceptable step.  Do this by fitting a model
% function through the existing data.  For the first pass this will be
% quadratic.  After that use a cubic function.  Lambda is the fraction of
% the Newton step that we should take.  (lambda = 1 means full Newton,
% while lambda = 0 means no step.)

% Set how many times to try before complaining.
jmax = 20;

% Mark the previous attempt at being the Newton step.
lambda1 = 1;
fL1 = fN;
drL = drN;

for j=1:1:jmax
    if(j == 1)
        % Use a quadratic model to project to the zero crossing.
        lambda = -dfdr/(2*(fL1-f-dfdr));
    else
        % Use a cubic model to better capture inflections.
        rhs1 = fL1 - f - lambda1*dfdr;
        rhs2 = fL2 - f - lambda2*dfdr;
        a = (rhs1/(lambda1*lambda1) - rhs2/(lambda2*lambda2))/(lambda1-lambda2);
        b = (-lambda2*rhs1/(lambda1*lambda1) + lambda1*rhs2/(lambda2*lambda2))/...
            (lambda1-lambda2);
        if(a == 0)
            lambda = -dfdr/b/2;
        else
            disc = b*b - 3*a*dfdr;
            if(disc < 0)
                lambda = 0.5*lambda1;
            else
                if(b <= 0)
                    lambda = (-b + sqrt(disc))/3/a;
                else
                    lambda = -dfdr/(b + sqrt(disc));
                end
            end
        end
    end
    
    % Don't accept lambda values that are too big or small.
    if(lambda > 0.5)
        lambda = 0.5;
    end
    if(lambda < 0.1)
        lambda = 0.1;
    end
    drL = lambda*drL;                   % Lambda step
    rL  = r + drL;                      % Lambda location
    fL  = P_crT(c,rL,T);         % Function at Lambda
    dfdrL = dPdr_crT(c,rL,T);    % Derivative at Lambda

    % See if the lambda step satisfies the decrease criterion.
    df_big_enough = ((f-fL) > -alpha*dfdr*drL);
    if(df_big_enough)
        % Accept the lambda step.
%         Iterations = j
        r = rL;
        dfdr = dfdrL;
        return
    end
    
    % Save values for next model.
    lambda2 = lambda1;
    lambda1 = lambda;
    fL2 = fL1;
    fL1 = fL;
end
if(j == jmax)
    disp('Backtrack not successfull in Liquid_Spinodal_cT')
end
