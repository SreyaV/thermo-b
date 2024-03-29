%% A1
clc;clear all;close all;
global N2 O2 Ar N P

N=3;
Setup_Air_Props;
P=100000;
%% Reboiler
% Set Ternary Air
x_out([N2, O2, Ar])=[0.025, 0.95, 0.025];

i=1;
for quality_out=0:0.1:1
    [T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, y_out, rho_f_out, rho_g_out, Q_req]=reboiler_bwd(x_out, quality_out);
    rho_g_out(i) = rho_g_out;
    x_data(i, :)=x_in;
    y_data(i, :)=y_out;
    T_data(i, 1:2)=[T_in, T_out];
    Q_data(i)=Q_req;
    quality_data(i)=quality_out;
    i=i+1;
end
%% Upward tray

i = 1;
for quality_out = 0:0.1:1
    i
    x_in = x_data(i,:);
    y_in = y_data(i,:);
    T_liq = T_data(i,1);
    T_vap = 89.1978;
    quality = quality_data(i);
    rg_in = rho_g_out(i);
    m_vap = quality;
    m_liq = 1;
    [x_out,y_out,T_liq_out,T_vap_out,m_vap_out,m_liq_in] = tray_upward(x_in,y_in,T_liq,T_vap,m_vap,m_liq);
    x_top(i,:) = x_out;
    y_top(i,:) = y_out;
    T_liq_sat(i,:) = T_liq_out;
    T_vap_sat(i,:) = T_vap_out;
    Quality(i) = m_vap_out/(m_liq_in+m_vap_out)
    i = i+1;

end

%% Condenser 

% Output reboiler- temperature, composition, density (liquid, vapour), complementary outlet phase, feed-specific heat transfer required 
clc;clear all;close all;
global N2 O2 Ar N P

N=3;
Setup_Air_Props;
P=100000;

% Set Ternary Air
y_feed([N2, O2, Ar])=[0.98, 0.01, 0.01];

i=1;
for quality_out=0:0.1:1
    [T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, x_out, rho_f_out, rho_g_out, Q_req]=condenser_bwd(y_feed, quality_out);
    rho_g_out(i) = rho_g_out;
    x_data(i, :)=x_out;
    y_data(i, :)=y_in;
    T_data(i, 1:2)=[T_in, T_out];
    Q_data(i)=Q_req;
    quality_data(i)=quality_out;
    i=i+1;
end

%% Downward tray

i = 1;
for quality_out = 0:0.1:1
    i
    x_in = x_data(i,:);
    y_in = y_data(i,:);
    T_vap = T_data(i,1);
    T_liq = 77.6402;
    quality = quality_data(i);
    m_vap = 1;
    m_liq = 1-quality;
    [x_out,y_out, T_liq_out,T_vap_out, m_vap_in, m_liq_out] = tray_downward(x_in, y_in, T_liq,T_vap,m_liq,m_vap);
    x_bot(i,:) = x_out;
    y_bot(i,:) = y_out;
    T_liq_sat(i,:) = T_liq_out;
    T_vap_sat(i,:) = T_vap_out;
    Quality(i) = m_vap_in/(m_liq_out+m_vap_in);
    i = i+1;

end

%% Upwards tray

clf
figure(1)
hold on 
color = ['g' 'b' 'r'];

%Below
plot(quality_data,x_top(:,1),'o-',color = 'g')
plot(quality_data,x_top(:,2),'o-',color = 'b')
plot(quality_data,10*x_top(:,3),'o-',color = 'r')

plot([0 1], [0 0], '-k')
plot([0 1], [0 0], '--k')
plot([0 1], [0 0], '+k')
plot([0 1], [0 0], 'ok')

plot(quality_data,y_top(:,1),'o--',color = 'g' )
plot(quality_data,y_top(:,2),'o--',color = 'b' )
plot(quality_data,10*y_top(:,3),'o--',color = 'r' )

%Above
plot(quality_data,x_data(:,1),'+-',color = 'g')
plot(quality_data,y_data(:,1),'+--',color = 'g' )
plot(quality_data,x_data(:,2),'+-',color = 'b')
plot(quality_data,y_data(:,2),'+--',color = 'b' )
plot(quality_data,10*x_data(:,3),'+-',color = 'r')
plot(quality_data,10*y_data(:,3),'+--',color = 'r' )

hold off
xlabel("Reboiler Outlet Quality (mass)")
ylabel("Mole Fraction")
title('Upwards Tray')
legend('Nitrogen ','Oxygen ','Argon (x10)','Liquid','Vapor','Below','Above');
plotfixer

%% 
%% Downwards tray

clf
figure(2)
hold on 
color = ['g' 'b' 'r'];

%Below
plot(quality_data,y_bot(:,1),'o--',color = 'g')
plot(quality_data,y_bot(:,2),'o--',color = 'b')
plot(quality_data,10*y_bot(:,3),'o--',color = 'r')

plot([0 1], [0 0], '-k')
plot([0 1], [0 0], '--k')
plot([0 1], [0 0], '+k')
plot([0 1], [0 0], 'ok')

plot(quality_data,x_bot(:,1),'o-',color = 'g' )
plot(quality_data,x_bot(:,2),'o-',color = 'b' )
plot(quality_data,10*x_bot(:,3),'o-',color = 'r' )

%Above
plot(quality_data,x_data(:,1),'+-',color = 'g')
plot(quality_data,y_data(:,1),'+--',color = 'g' )
plot(quality_data,x_data(:,2),'+-',color = 'b')
plot(quality_data,y_data(:,2),'+--',color = 'b' )
plot(quality_data,10*x_data(:,3),'+-',color = 'r')
plot(quality_data,10*y_data(:,3),'+--',color = 'r' )

hold off
xlabel("Condenser Outlet Quality (mass)")
ylabel("Mole Fraction")
title('Downwards Tray')
legend('Nitrogen ','Oxygen ','Argon (x10)','Liquid','Vapor','Below','Above');
plotfixer



%% Upwards tray Temperatures
clf
figure(2)
plot(quality_data(1:end),T_liq_sat,"o-",color='blue')
hold on
plot(quality_data(1:end),T_vap_sat,"o-",color='black')
plot(quality_data,T_data(:,2),'+--',color='blue')
legend('Liquid Above', 'Outlet Liquid & Vapor','Vapor Below')
ylabel('Temperature (K)')
xlabel('Reboiler Outlet Quality (mass)')
plotfixer

%% 

%% Upwards tray

clf
figure(1)
hold on 
color = ['g' 'b' 'r'];

%Below
plot(quality_data,x_top(:,1),'o-',color = 'g')
plot(quality_data,x_top(:,2),'o-',color = 'b')
plot(quality_data,10*x_top(:,3),'o-',color = 'r')

plot([0 1], [0 0], '-k')
plot([0 1], [0 0], '--k')
plot([0 1], [0 0], '+k')
plot([0 1], [0 0], 'ok')

plot(quality_data,y_top(:,1),'o--',color = 'g' )
plot(quality_data,y_top(:,2),'o--',color = 'b' )
plot(quality_data,10*y_top(:,3),'o--',color = 'r' )

%Above
plot(quality_data,x_data(:,1),'+-',color = 'g')
plot(quality_data,y_data(:,1),'+--',color = 'g' )
plot(quality_data,x_data(:,2),'+-',color = 'b')
plot(quality_data,y_data(:,2),'+--',color = 'b' )
plot(quality_data,10*x_data(:,3),'+-',color = 'r')
plot(quality_data,10*y_data(:,3),'+--',color = 'r' )

hold off
xlabel("Reboiler Outlet Quality (mass)")
ylabel("Mole Fraction")
title('Upwards Tray')
legend('Nitrogen ','Oxygen ','Argon (x10)','Liquid','Vapor','Below','Above');
plotfixer

%% Downwards tray Temperatures
clf
figure(2)
hold on
plot(quality_data(1:end),T_vap_sat,"+--",color='blue')
plot(quality_data(1:end),T_liq_sat,"o-",color='black')
plot(quality_data,T_data(:,2),'o-',color='blue')
legend('Vapor Below', 'Outlet Liquid & Vapor','Liquid Above')
ylabel('Temperature (K)')
xlabel('Condenser Outlet Quality (mass)')
plotfixer


%% Upwards Tray outlet quality
clf
plot(quality_data,Quality,'o-r')
ylabel('Tray Outlet Quality (mass)')
xlabel('Reboiler Outlet Quality (mass)')
plotfixer

%% Downwards Tray outlet quality
clf
plot(quality_data,Quality,'o-r')
ylabel('Tray Outlet Quality (mass)')
xlabel('Condenser Outlet Quality (mass)')
plotfixer