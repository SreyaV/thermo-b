function [T rl rv y] = Fast_Bubble_cP(x,P,varargin)
% Interpolation wrapper for Bubble_cP that is faster for binary mixtures.
% For binary mixtures, and a given pressure, this records all the past
% results and starts to build interpolation arrays. It tries getting a
% result by interpolating first. If the interpolated result is close
% enough, it returns that. Otherwise, it calls Bubble_cP.
% Original function by Kevin Shan, ME370B class of 2010.
% Rewritten for general use by C.F. Edwards, 2-26-12

% Save things we need between calls.
global nShortcutted nInArray Array
% Get the tolerance info from outside.
global toler O2 N2

% Check the number of species in the mixture
if length(x) ~= 2
    % Doesn't currently support non-binary mixtures.
    if nargin == 6
        [T rl rv y] = Bubble_cP(x,P,varargin);
    else
        [T rl rv y] = Bubble_cP(x,P);
    end        
    return
end

% Check to see if this is the first call.
if isempty(Array)
    % First time the routine has been called.
    % Initialize the array and number of pressures it contains.
    nPressures = 0;
    Array(1).P = 0;
else
    nPressures = length(Array);
end

% Test to see if we have this pressure.
for i=1:1:length(Array)
    if P == Array(i).P
        isInArray = true;
        index = i;
        break
    else
        isInArray = false;
    end
end

% If the pressure is not in the array, we need to add it.
% Note that the intent is to handle a few fixed pressures, not a continuous
% range of pressure.  If you want that, you should re-write this to do a
% two-dimensional interpolation in pressure and composition.
if ~isInArray
    % Increment the index for a new pressure value.
    index = nPressures + 1;
    Array(index).P = P;
    % Make coarse interp arrays to start
    fprintf('Fast_Bubble_cP: New pressure P = %0.1f bar\n', P/1e5);
    nPts = 5;
    Array(index).N2l = linspace(0.05,1-.05,nPts)';
    % Add in the two end points.
    Array(index).N2l = [0; Array(index).N2l; 1];
    nPts = nPts+2;
    Array(index).T   = zeros(nPts,1);
    Array(index).rl  = zeros(nPts,1);
    Array(index).rv  = zeros(nPts,1);
    Array(index).N2v = zeros(nPts,1);
    % Populate the interpolation arrays
    fprintf('Init interp arrays. N2 = ');
    % Do the pure O2 point.
    ii = 1;
    N2l = Array(index).N2l(ii);
    fprintf('%0.2f, ',N2l);
    x_test = [N2l 1-N2l];
    % Calculate
    [T rl rv] = Saturation_iP(O2,P);
    y = x_test;
    % Store
    Array(index).T(ii)   = T;
    Array(index).rl(ii)  = rl;
    Array(index).rv(ii)  = rv;
    Array(index).N2v(ii) = y(1);
    for ii = 2:nPts-1
        % Setup
        N2l = Array(index).N2l(ii);
        fprintf('%0.2f, ',N2l);
        x_test = [N2l 1-N2l];
        % Calculate
        [T rl rv y] = Bubble_cP(x_test,P);
        % Store
        Array(index).T(ii)   = T;
        Array(index).rl(ii)  = rl;
        Array(index).rv(ii)  = rv;
        Array(index).N2v(ii) = y(1);
    end
    % Do the pure N2 point.
    ii = ii+1;
    N2l = Array(index).N2l(ii);
    fprintf('%0.2f, ',N2l);
    x_test = [N2l 1-N2l];
    % Calculate
    [T rl rv] = Saturation_iP(N2,P);
    y = x_test;
    % Store
    Array(index).T(ii)   = T;
    Array(index).rl(ii)  = rl;
    Array(index).rv(ii)  = rv;
    Array(index).N2v(ii) = y(1);
    fprintf('\nFast_Bubble_cP: interpolation arrays initialized\n');
    nShortcutted(index) = 0;
    nInArray(index) = length(Array(index).N2l);
end

% So the pressure is now in our storage array and index points to it.
% Interpolate to get a guess first.
% Get N2 molefraction from the input composition
N2l = x(1);
% Clip if out of bounds.
if N2l > 1
    N2l = 1;
elseif N2l < 0
    N2l = 0;
end
x = [N2l  1-N2l];
% Interpolate output
T   = interp1(Array(index).N2l, Array(index).T,   N2l, 'spline'); % Temperature
rl  = interp1(Array(index).N2l, Array(index).rl,  N2l, 'spline'); % Liquid density
rv  = interp1(Array(index).N2l, Array(index).rv,  N2l, 'spline'); % Vapor density
N2v = interp1(Array(index).N2l, Array(index).N2v, N2l, 'spline'); % Vapor N2 molefraction
% Clip if out of bounds.
if N2v > 1
    N2v = 1;
elseif N2v < 0
    N2v = 0;
end
y  = [N2v  1-N2v];
% Check pressure and chem. potential tolerances
P_interp =   [P_crT(x,rl,T)       P_crT(y,rv,T)];
mul_interp = [mui_icrT(1,x,rl,T)  mui_icrT(2,x,rl,T)];
muv_interp = [mui_icrT(1,y,rv,T)  mui_icrT(2,y,rv,T)];
% Tolerance
TOLER = 10*toler;
dT = TOLER*300;
dmu = abs(muv_interp(1)-mui_icrT(1,y,rv,T+dT));
if (abs((P_interp(1)-P)/P) < TOLER)...
        && (abs((P_interp(2)-P)/P) < TOLER)...
        && (abs(mul_interp(1)-muv_interp(1)) < dmu)...
        && (abs(mul_interp(2)-muv_interp(2)) < dmu)
    % Within tolerance.
    nShortcutted(index) = nShortcutted(index) + 1;
    return;
end

% Interpolation was not successful.  Call Bubble_cP
try
    % Give Bubble_cP a starting guess based on our interpolation
    [T rl rv y] = Bubble_cP(x,P,T,y,rl,rv);
    % If that didn't work, then it might return a nonsensical result.
    % In that case, try Bubble_cP without a starting guess
    if (min(y) < 0) || ~isreal(T)
        [T rl rv y] = Bubble_cP(x,P);
    end
catch ME
    % This is in case Bubble_cP exited abnormally when we gave it a
    % starting guess
    if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
        [T rl rv y] = Bubble_cP(x,P);
    else
        rethrow(ME)
    end
end

% Add this to our interp arrays
N2v = y(1);
% Check if this is already in our interp arrays. This happens when the
% result returned by Bubble_cP is out of tolerance, so that future
% interpolations based on that result are also out of tolerance (even if a
% future call is at exactly the same point).
if find(Array(index).N2l == N2l)
    % Do nothing...assume the older value is correct
    return
end
% Add this result to our interp arrays and sort
[Array(index).N2l, sortIndex] = sort([Array(index).N2l; N2l]);
Array(index).T   = [  Array(index).T; T];
Array(index).T   = Array(index).T(sortIndex);
Array(index).rl  = [ Array(index).rl; rl];
Array(index).rl  = Array(index).rl(sortIndex);
Array(index).rv  = [ Array(index).rv; rv];
Array(index).rv  = Array(index).rv(sortIndex);
Array(index).N2v = [Array(index).N2v; N2v];
Array(index).N2v = Array(index).N2v(sortIndex);
% Update the length value for return to calling script.
nInArray(index) = length(Array(index).N2l);
