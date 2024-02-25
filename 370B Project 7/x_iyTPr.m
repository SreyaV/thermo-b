function x = x_iyTPr(isp,y,T,P,rl,rv,xstart)
% Return the liquid molefractions that match the chemical potentials of the
% vapor for all but the highest molefraction species.  Set that one (isp)
% to be the remainder of the mixture.  (It will be tested in an outer loop.)
% C.F. Edwards, 2-16-10

global toler

% Assign some storage used inside loops (to avoid dynamic allocation).
muvi  = zeros(length(y),1);
intol = zeros(length(y),1);

% Use a supplied starting point.
x = xstart;

% Find the target chemical potentials (set by the vapor).
rv = rv_cTP(y,T,P,rv);
for i=1:1:length(y)
    if(i ~= isp)
        muvi(i) = mui_icrT(i,y,rv,T);
    else
        muvi(i) = 0;
    end
end

jmax = 50;
for j=1:1:jmax
    for i=1:1:length(x)
        if(i ~= isp)
            % Use N-R to get the ith liquid molefraction to match the chemical
            % potential of the vapor while holding the others constant and setting
            % isp as the remainder.  
            
            % Normalize by adjusting the isp molefraction.
            x(isp) = 0;
            x(isp) = 1 - sum(x);

            % Get the density and chemical potential for this composition.
            rl = rl_cTP(x,T,P,rl);
            mul = mui_icrT(i,x,rl,T);
            
            % Find the deviation and use NR to update the species.
            f = mul - muvi(i);                     
            dmudc = dmuvdci_icTPrmu(i,x,T,P,rl,mul);
            x(i) = x(i) - f/dmudc;
            if(x(i) < 0)
                x(i) = 0;
            end
            if(x(i) > 1)
                x(i) = 1;
            end
        end
    end
    x(isp) = 0;
    x(isp) = 1 - sum(x);

    % Test this composition for the correct joint molefractions.
    for i=1:1:length(x)
        if(i ~= isp)
            rl = rl_cTP(x,T,P,rl);
            mul = mui_icrT(i,x,rl,T);
            f = mul - muvi(i);
            intol(i) = (abs(f/muvi(i)) < toler);
        else
            intol(i) = 1;
        end
    end
    if(sum(intol) == length(y))
        return
    end
end
disp('Fell off the end of NR loop in x_iyTPr')
