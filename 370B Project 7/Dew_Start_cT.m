function [P rg rf x] = New_Dew_Start_cT(y,T)
% Return a set of starting values for a dew point calculation.  Use
% a direct ideal solution approach where possible (low temp) and an
% interpolation approach where that fails.  The first time you do the
% interpolation, save it in non-volatile storage for possible future use.
% C.F. Edwards, 2-20-12

global N2 O2 Ar
persistent stored

% Copy the composition to avoid overwriting.
c = y;

N = length(c);  % Find out how many components.
if ~((N == 2)||(N == 3))
    disp('Number of components must be 2 or 3 in Dew_Start_cT')
    return
end

% Try the direct approach.
Tmin = 60;
Tmax = 120;
if T <= Tmax
    [P rg rf x] = Ideal_Dew_cT(c,T);
    % Successful?
    if(P ~= 0)
        return
    end
end

% Must construct a value by interpolation.  Will want to store the result
% in case we ask for another starting point with the same composition.
c_toler = 1e-2;         % How close do we have to be on composition?
in_range = zeros(N,1);  % Logic flag array for components.
if ~isempty(stored)
    % Some composition(s) has (have) been done already.
    % Find out if we are close enough to a stored composition.
    for j=1:1:stored.n
        % Check each independent component.
        for i=1:1:N
            in_range(i) = ...
                ((1-c_toler)*stored.c(i,j) < c(i))&&(c(i) < (1+c_toler)*stored.c(i,j));
        end
        % Assemble a logic flag and test it.
        if sum(in_range)==N
            % We have this composition stored.  Use it and return.
            % Interpolate to find the data at the temperature specified.
            rf = interp1(stored.Tspan(:,j),stored.rfspan(:,j),T,'spline');
            rg = interp1(stored.Tspan(:,j),stored.rgspan(:,j),T,'spline');
            P  = interp1(stored.Tspan(:,j),stored.Pspan(:,j),T,'spline');
            x(N2) = interp1(stored.Tspan(:,j),stored.N2span(:,j),T,'spline');
            x(O2) = interp1(stored.Tspan(:,j),stored.O2span(:,j),T,'spline');
            if(N == 3)
                x(Ar) = interp1(stored.Tspan(:,j),stored.Arspan(:,j),T,'spline');
            end
            x = real(x/sum(x));
            return
        end
    end
    % If we got to here, no entry matches.
    % Make a new entry and then go and solve it.
    stored.n = stored.n+1;
    stored.c(:,stored.n) = c;
else
    % If we got to here, there are no prior attempts.  Define the variables
    % and then go and solve it.
    stored.n = 1;
    stored.c(:,1) = c;
end

% Get the P-rho inflection point to see roughly where the top of the dome
% should be.
[T_Pri r_Pri] = Pr_Inflection_c(c);
P_Pri = P_crT(c,r_Pri,T_Pri);

% Step through collecting estimates of the dew point until we can go no
% further.  Limit the range of attempts to less than Tmax.
dT = 5;
i = 1;
% Use a lower temperature limit for the fit since otherwise the answer
% pulls in too much on the liquid density.
Tmaxfit = 100;
for TT=Tmin:dT:Tmaxfit
    TTd = TT;
    [Pdew rgdew rfdew xdew] = Ideal_Dew_cT(c,TT);
    if(Pdew == 0)
        break
    end
    % Save the values.
    Tdewi(i)   = TT;
    Pdewi(i)   = Pdew;
    rfdewi(i)  = rfdew;
    rgdewi(i)  = rgdew;
    xdewi(i,:) = xdew;
    i = i+1;
end

% Merge these into a set of arrays that span the density range.
% Use the inflection point--shifted as needed--to set the location to which
% the interpolation functions will point.
rfspan = [rfdewi 1.51*r_Pri];
rgspan = [rgdewi 0.71*r_Pri];
Tspan = [Tdewi  T_Pri];
Pspan = [Pdewi  0.987*P_Pri];
N2span = [xdewi(:,N2)' 0.91*c(N2)];
O2span = [xdewi(:,O2)' (1/.91)*c(O2)];
if(N == 3)
    Arspan = [xdewi(:,Ar)' c(Ar)];
end

% Store these for fast response next time this composition is chosen.
stored.rfspan(:,stored.n) = rfspan';
stored.rgspan(:,stored.n) = rgspan';
stored.Tspan(:,stored.n) = Tspan';
stored.Pspan(:,stored.n) = Pspan';
stored.N2span(:,stored.n) = N2span';
stored.O2span(:,stored.n) = O2span';
if(N == 3)
    stored.Arspan(:,stored.n) = Arspan';
end

% Interpolate to find the data at the temperature specified.
rf = interp1(Tspan,rfspan,T,'spline');
rg = interp1(Tspan,rgspan,T,'spline');
P  = interp1(Tspan,Pspan,T,'spline');
x(N2) = interp1(Tspan,N2span,T,'spline');
x(O2) = interp1(Tspan,O2span,T,'spline');
if(N == 3)
    x(Ar) = interp1(Tspan,Arspan,T,'spline');
end    
x = real(x/sum(x));

% % The remainder is for debugging.  Uncomment to see plots.
% 
% % Interplote the data across the range.
% Tint  = linspace(min(Tspan),max(Tspan),50);
% rfint = interp1(Tspan,rfspan,Tint,'spline');
% rgint = interp1(Tspan,rgspan,Tint,'spline');
% Pint  = interp1(Tspan,Pspan,Tint,'spline');
% N2int = interp1(Tspan,N2span,Tint,'spline');
% O2int = interp1(Tspan,O2span,Tint,'spline');
% if(N == 3)
%     Arint = interp1(Tspan,Arspan,Tint,'spline');
% end
% 
% figure(1)
% clf
% hold on
% plot(rgspan,Tspan,'ko')
% plot(rgint,Tint,'k-')
% plot(r_Pri,T_Pri,'ro')
% plot(rg,T,'kd')
% hold off
% xlabel('Density (kg/m3)')
% ylabel('Temperature (K)')
% plotfixer
% 
% figure(2)
% clf
% hold on
% plot(rgspan,Pspan/1e6,'ko')
% plot(rgint,Pint/1e6,'k-')
% plot(r_Pri,P_Pri/1e6,'ro')
% plot(rg,P/1e6,'kd')
% hold off
% xlabel('Density (kg/m3)')
% ylabel('Pressure (MPa)')
% plotfixer
% 
% figure(3)
% clf
% hold on
% plot(rgspan,N2span,'o','Color',[0 .5 0])
% plot(rgspan,O2span,'bo')
% plot(rgspan,10*Arspan,'ro')
% plot(rgint,N2int,'-','Color',[0 .5 0])
% plot(r_Pri,c(N2),'kx')
% plot(rgint,O2int,'b-')
% plot(r_Pri,c(O2),'kx')
% plot(rgint,10*Arint,'r-')
% plot(r_Pri,10*c(Ar),'kx')
% plot(rg,x(N2),'d','Color',[0 .5 0])
% plot(rg,x(O2),'bd')
% plot(rg,10*x(Ar),'rd')
% hold off
% xlabel('Density (kg/m3)')
% ylabel('Mole Fraction')
% legend('N2','O2','Ar',0)
% plotfixer
% 
% figure(4)
% clf
% hold on
% plot(rgspan,rfspan,'ko')
% plot(rgint,rfint,'k-')
% plot(r_Pri,r_Pri,'ro')
% plot(rg,rf,'kd')
% hold off
% xlabel('Vapor Density (kg/m3)')
% ylabel('Liquid Density (kg/m3)')
% plotfixer
