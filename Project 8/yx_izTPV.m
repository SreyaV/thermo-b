function [y x] = yx_izTPV(isp,z,T,P,V,ystart,rlstart,rvstart)
% Return the vapor and liquid molefractions that match in chemical potential 
% for all but the highest molefraction species.  Set that one (isp)
% to be the remainder of the mixture.  (It will be tested in an outer
% loop.)  Use the specified molar quality V to couple the species changes. 
% C.F. Edwards, 2-16-10

global toler

% Assign some storage used inside loops (to avoid dynamic allocation).
x = zeros(1,length(z));
intol = zeros(length(z),1);

% The x increment is coupled to yinc by the molar quality because of the 
% constraint on the total amount of species i.  The species is transferred
% from one phase to the other, hence the increment goes in opposite
% directions.
yinc = 1e-4;            % Use a small increment for y derivatives.
xinc = -yinc*V/(1-V);   % Set the x increment to enforce species continuity.

% Set the starting vapor composition.
y = ystart;             

% Make a liquid that obeys the y,z,V constraints.
for i=1:1:length(z)
    x(i) = (z(i) - y(i)*V)/(1-V);    
end

% If the liquid has negative mole numbers, make them small and adjust vapor
% to fit the constraints.
for i=1:1:length(z)
    if(i ~= isp)
        if(x(i) <= 0)
            % Set to zero.
            x(i) = 0;
        end
    else
        % Clear the isp value.
        x(i) = 0;
    end
end
% Use isp for closure.
x(isp) = 1 - sum(x);
% Make the vapor obey the z,V constraints.
for i=1:1:length(z)
    y(i) = (z(i) - x(i)*(1-V))/V;    
end

rv = rvstart;
rl = rlstart;

jmax = 50;
for j=1:1:jmax
    for i=1:1:length(z)
        if(i ~= isp)
            % Use N-R to get the ith vapor molefraction to match the chemical
            % potential of the liquid while holding the others constant and setting
            % isp as the remainder.
            rv = rv_cTP(y,T,P,rv);
            muiv = mui_icrT(i,y,rv,T);
            rl = rl_cTP(x,T,P);
            muil = mui_icrT(i,x,rl,T);  
            f = muiv - muil;
            
            % Use NR to adjust the composition to match the chemical potentials.
            % Use a forward difference for the numerical derivative.
            yihigh     = y(i) + yinc;
            yhigh      = y;
            yhigh(i)   = yihigh;
            yhigh(isp) = 0;
            yhigh(isp) = 1 - sum(yhigh);
            rv = rv_cTP(yhigh,T,P,rv);
            muivhigh = mui_icrT(i,yhigh,rv,T);
            xihigh     = x(i) + xinc;
            xhigh      = x;
            xhigh(i)   = xihigh;
            xhigh(isp) = 0;
            xhigh(isp) = 1 - sum(xhigh);
            rl = rl_cTP(xhigh,T,P,rl);
            muilhigh = mui_icrT(i,xhigh,rl,T);
            
            dfdyi = ((muivhigh-muilhigh)-(muiv-muil))/(yinc);

            ylast = y;
            xlast = x;
            dyi = -f/dfdyi;
            dxi = -dyi*V/(1-V);
            
            y(i)   = y(i) + dyi;
            if(y(i) < 0)
                y(i) = 0;
            end
            if(y(i) > 1)
                y(i) = 1;
            end
            y(isp) = 0;
            y(isp) = 1 - sum(y);

            x(i)   = x(i) + dxi;
            if(x(i) < 0)
                x(i) = 0;
            end
            if(x(i) > 1)
                x(i) = 1;
            end
            x(isp) = 0;
            x(isp) = 1 - sum(x);
        end
    end

    % Test this composition for the correct joint molefractions.
    for i=1:1:length(z)
        if(i ~= isp)
            rv = rv_cTP(y,T,P,rv);
            muiv = mui_icrT(i,y,rv,T);
            rl = rl_cTP(x,T,P,rl);
            muil = mui_icrT(i,x,rl,T);
            f = muiv - muil;
            intol(i) = ((abs(f/muil) < toler)&&...
                ((abs(y(i)-ylast(i))/y(i)) < toler)&&...
                ((abs(x(i)-xlast(i))/x(i)) < toler) );
        else
            intol(i) = 1;
        end
    end
    if(sum(intol) == length(z))
        break
    end
end

if(j == jmax)
    z
    T
    P
    disp('Fell off the end of NR loop in yx_izTPV')
    error('Fell off the end of NR loop in yx_izTPV')
end
y = real(y);
x = real(x);

