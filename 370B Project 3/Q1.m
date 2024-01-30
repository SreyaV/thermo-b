%% Close Everything
clc;clear;close all;
fprintf('\nClosing everything...\n')
%% Declare variables
fprintf('\nDeclaring Variables...\n')
T=[1600,1800,2000];                         % Test cases given in question 1
PRlim=[100 150 200];                        % PR Range inferred from Chris' lecture PPT for corresponding Temperatures
Reso=2.5;                                     % Set resolution for lines on graph (must be factor of all PRLim)
%% Wrapper loop
for temp=1:1:length(T)
    fprintf('\nCalculating for %dK', T(temp))                 % Disp which temperature is being used for calculation
    PR=10:Reso:PRlim(temp);
    for j=1:1:length(PR)
        [eff(j, temp), w(j, temp)]=gas_turbine(PR(j),T(temp));
        fprintf('.')
    end
end
%% Plot
fprintf('\nPlotting...\n');
figure()
hold on
for i=1:1:length(T)
    for j=1:1:length(eff)
        if eff(j, i)~=0
            x(j)=eff(j, i)*100;
            y(j)=w(j, i)/1000;
        end
    end
    plot(y, x);
    clear x y;
end
legend('1600K','1800K','2000K')
plotfixer
fprintf('\nDone!\n')