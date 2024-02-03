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
figure(2);
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



plot(w(1, 1)/1000, eff(1, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(1, 1)/1000, eff(1, 1)*100, ' 10', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(w(5, 1)/1000, eff(5, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(5, 1)/1000, eff(5, 1)*100, '20', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(w(13, 1)/1000, eff(13, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(13, 1)/1000, eff(13, 1)*100, ' 40', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(w(21, 1)/1000, eff(21, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(21, 1)/1000, eff(21, 1)*100, ' 60', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(w(29, 1)/1000, eff(29, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(29, 1)/1000, eff(29, 1)*100, ' 80', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
plot(w(37, 1)/1000, eff(37, 1)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty 
text(w(37, 1)/1000, eff(37, 1)*100, ' 100', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%
plot(w(1, 2)/1000, eff(1, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(1, 2)/1000, eff(1, 2)*100, ' 10', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(5, 2)/1000, eff(5, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(5, 2)/1000, eff(5, 2)*100, ' 20', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(13, 2)/1000, eff(13, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(13, 2)/1000, eff(13, 2)*100, ' 40', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(23, 2)/1000, eff(23, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(23, 2)/1000, eff(23, 2)*100, ' 65', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(37, 2)/1000, eff(37, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(37, 2)/1000, eff(37, 2)*100, ' 100', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(47, 2)/1000, eff(47, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(47, 2)/1000, eff(47, 2)*100, ' 125', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(57, 2)/1000, eff(57, 2)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Large empty circle
text(w(57, 2)/1000, eff(57, 2)*100, ' 150', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%

plot(w(1, 3)/1000, eff(1, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(1, 3)/1000, eff(1, 3)*100, ' 10', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(5, 3)/1000, eff(5, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(5, 3)/1000, eff(5, 3)*100, ' 20', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(13, 3)/1000, eff(13, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(13, 3)/1000, eff(13, 3)*100, ' 40', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(23, 3)/1000, eff(23, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(23, 3)/1000, eff(23, 3)*100, ' 65', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(37, 3)/1000, eff(37, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(37, 3)/1000, eff(37, 3)*100, ' 100', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(57, 3)/1000, eff(57, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(57, 3)/1000, eff(57, 3)*100, ' 150', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot(w(77, 3)/1000, eff(77, 3)*100, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerFaceColor', 'none');
text(w(77, 3)/1000, eff(77, 3)*100, ' 200', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlabel('Air-Specific Work, ')
legend('1600K','1800K','2000K')
hold off;
% % Coordinates of the point
% px = 5;
% py = 10;
% 
% % Plot the point as a large empty circle
% figure; % Opens a new figure window
% plot(x, y, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'); % Large empty circle
% 
% % Keep the plotq4_exergy_plots.m open to add a label
% hold on;
% 
% % Add a label to the point
% text(x, y, ' Your Label', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% 

plotfixer
fprintf('\nDone!\n')



