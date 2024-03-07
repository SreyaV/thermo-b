% Generate a 3-D T-s-x surface plot for two component N2 and O2.
% Load a pre-computed and stored data space for speed.
% Use Matlab's inherent graphics to rotate it around to view the space.
% C.F. Edwards, 2-14-10

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

format compact
fprintf('\n************************************************************\n')

% Load a file to save recalculating things up to here.
load Tsx_Data

% We will build the surface as we go.  Start a figure and put it on hold.
figure(1)
clf
hold on
ylabel('Nitrogen Mole Fraction','rotation',0)
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
zlabel('Temperature (K)')
view([-10 40])
axis([2 7 0 1 10 300])
grid on
drawnow

% Tell the user the composition planes to be shown.
xN2_planes = clist

% Show the P-rho critical locus on the plot.
plot3(sc/1e3,xN2,Tc,'kx','LineWidth',2)

% Add domes at each plane.
for i=1:1:NCS
    plot3(sdome(i,:)/1e3,xN2(i)*ones(1,length(sdome)),Tdome(i,:),'k-','LineWidth',2)
    drawnow
end

% Add dew and bubble points.
for i=1:1:NCS
    plot3(sdew(i,:)/1e3,xN2(i)*ones(1,NPSSC),Tdew(i,:),'b.','LineWidth',2)
    plot3(sbub(i,:)/1e3,xN2(i)*ones(1,NPSSC),Tbub(i,:),'b.','LineWidth',2)
    drawnow
end

% Put splines through the bubble and dew data in the composition direction.
[row col] = size(Tdew);
for j=1%:1:col
    plot3(sdewSpline(j,:)/1e3,xN2Spline,TdewSpline(j,:),'b-','LineWidth',2)
    plot3(sbubSpline(j,:)/1e3,xN2Spline,TbubSpline(j,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist

% Add isobars.
for i=1:1:NCS
    % Do the points below the dome.
    for j=1%:1:col
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
    % Now do the pressures that do not have dew points (above dome).
    for j=j+1%:1:NPS
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
end

%% 10 trays
global N2 O2 Ar N P
N_trays = 10;
reboiler_quality = 0.7;
N = 2;
N2 = 1;
O2 = 2;
P = 1e5;
Setup_Air_Props;
[feed,tray,out,x_out] = stripper_wrap(reboiler_quality,N_trays);

%% Plot flow rate (N_trays = 10)
clf
vapor_flow = [out.vap out.tot];
liquid_flow = [out.liq out.tot];
tot_flow = [out.tot];
for k=1:N_trays
    vapor_flow = [vapor_flow tray(k).vap tray(k).tot];
    liquid_flow = [liquid_flow tray(k).liq tray(k).tot ];
    tot_flow = [tot_flow tray(k).tot];
end
vapor_flow = [vapor_flow feed.vap feed.vap+1.8];
liquid_flow = [liquid_flow feed.liq feed.liq];
tot_flow = [tot_flow feed.tot+1.8];

plot(vapor_flow,"o--g")
hold on
plot(liquid_flow,"o--b")
plot([2 4 6 8 10 12 14 16 18 20 22 24],tot_flow,'ok')

legend("Vapor Flow","Liquid Flow","Total Flow")
xlabel("Stage")
xticks([2 4 6 8 10 12 14 16 18 20 22 24])
xticklabels({'R','1','2','3','4','5','6',"7",'8','9','10','F'})
ylabel("Semi-Extensive Flow Rate (mol/mol of product)")
title("10 stages, 0.7 Reboiler Quality ")
hold off
plotfixer

%% Plot T-s-x

% Trays
hold on
for j=1:N_trays
    plot3([tray(j).xS tray(j).N2s tray(j).yS]/1e3,[tray(j).x tray(j).N2 tray(j).y],[tray(j).T tray(j).T tray(j).T],'x-g','LineWidth',2)
    drawnow
end
for j=1:N_trays-1
    plot3([tray(j).N2s tray(j+1).N2s]/1e3,[tray(j).N2 tray(j+1).N2],[tray(j).T tray(j+1).T],'o--r','LineWidth',2)
    drawnow
end

% Out
plot3([out.xS out.N2s out.yS]/1e3,[out.x out.N2 out.y],[out.T out.T out.T],'x-g','LineWidth',2)
plot3([tray(1).N2s out.N2s]/1e3,[tray(1).N2 out.N2],[tray(1).T out.T],'o--r','LineWidth',2)

% Feed
plot3([feed.xS feed.N2s feed.yS]/1e3,[feed.x feed.N2 feed.y],[feed.T feed.T feed.T],'x-g','LineWidth',2)
plot3([feed.N2s tray(N_trays).N2s]/1e3,[feed.N2 tray(N_trays).N2],[feed.T tray(N_trays).T],'o--r','LineWidth',2)

%Flash
P = 1e5;
z_feed=[0.71 0.29];
s1 = s_crT(z_feed,rv_cTP(z_feed,300,P),300);
s2 = s_crT(z_feed,rv_cTP(z_feed,300,P*200),300);
s3 = s_crT(z_feed,rl_cTP(z_feed,130,P*200),130);
plot3([s1 s2 s3 feed.N2s]/1e3,[z_feed(1) z_feed(1) z_feed(1) feed.N2],[300 300 130 feed.T],'o--r','LineWidth',2)
%plot3([3 3],[feed.N2 feed.N2],[feed.T 150],'--r','LineWidth',2)

title("10 stages, 0.7 reboiler quality")

% Finished with the dynamic plot.  Close it out.
hold off

%%
a= [[0 -1 0 0.5];[0 0 1 -0.5];[1 0 0 8];[0 0 0 1]];
