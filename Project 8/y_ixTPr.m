function y = y_ixTPr(isp,x,T,P,rl,rv,ystart)
% Return the vapor molefractions that match the chemical potentials of the
% liquid for all but the highest molefraction species.  Set that one (isp)
% to be the remainder of the mixture.  (It will be tested in an outer loop.)
% C.F. Edwards, 2-21-10

global toler

% Assign some storage used inside loops (to avoid dynamic allocation).
muli  = zeros(length(x),1);
muv   = zeros(length(x),1);
f     = zeros(length(x),1);
intol = zeros(length(x),1);

% Use the supplied starting point.
y = Clean_c(ystart);
if sum(y) == 0
    disp('Invalid composition given to y_ixTPr')
    return
end

% Find the target chemical potentials (set by the liquid).
rv = rv*ones(3);
for i=1:1:length(x)
    if(i ~= isp)
        muli(i) = mui_icrT(i,x,rl,T);
        rv(i) = rv_cTP(y,T,P,rv(i));
        muv(i) = mui_icrT(i,y,rv(i),T);
        f(i) = muv(i) - muli(i);
    else
        muli(i) = 0;
    end
end

jmax = 50;
for j=1:1:jmax
    for i=1:1:length(y)
        if(i ~= isp)
            % Use N-R to get the ith vapor molefraction to match the chemical
            % potential of the liquid while holding the others constant and setting
            % isp as the remainder.
            ylast = y(i);
            
            dmudc = dmuvdci_icTPrmu(i,y,T,P,rv(i),muv(i));
            dy = -f(i)/dmudc;
            limit = 0.1;
            if(abs(dy) > limit)
                dy = sign(dy)*limit;
            end
            y(i) = y(i) + dy;
            % Bisect if N-R runs out of bounds.
            if(y(i) < 0)
                y(i) = (ylast+0)/2;
            end
            if(y(i) > 1)
                y(i) = (1+ylast)/2;
            end
            y(isp) = 0;
            y(isp) = 1 - sum(y);
            y = real(y);
            if(y(isp) < 0)
                y(isp) = toler;
                y(i) = 0;
                y(i) = 1 - sum(y);
            end
            if(y(isp) > 1)
                y(isp) = 1-toler;
                y(i) = 0;
                y(i) = 1 - sum(y);
            end
            y = real(y);
        end
    end

    y(isp) = 0;
    y(isp) = 1 - sum(y);
    y = real(y);
    
    % Test this composition for the correct joint molefractions.
    for i=1:1:length(y)
        if(i ~= isp)
            rv(i) = rv_cTP(y,T,P,rv(i));
            muv(i) = mui_icrT(i,y,rv(i),T);
            f(i) = muv(i) - muli(i);
            intol(i) = (abs(f(i)/muli(i)) < toler);
        else
            intol(i) = 1;
        end
    end
    if(sum(intol) == length(x))
        return
    end
end
j
x
disp('Fell off the end of NR loop in y_ixTPr')
