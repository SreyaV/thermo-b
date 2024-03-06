%% Setup Reboiler

% Output reboiler- temperature, composition, density (liquid, vapour), complementary outlet phase, feed-specific heat transfer required 
clc;clear all;close all;
global N2 O2 Ar N P

N=3;
Setup_Air_Props;
P=100000;

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

i = 9;
for quality_out = 0.8%0:0.1:1
    i
    x_in = x_data(i,:);
    y_in = y_data(i,:);
    T_liq = T_data(i,1);
    T_vap = 89.1978;
    quality = quality_data(i);
    rg_in = rho_g_out(i);
    m_vap = quality;
    m_liq = 1-quality;
    [x_out,y_out,T_liq_out,T_vap_out] = tray_upward(x_in,y_in,T_liq,T_vap,m_vap,m_liq);
    x_top(i,:) = x_out;
    y_top(i,:) = y_out;
    T_liq_sat(i,:) = T_liq_out;
    T_vap_sat(i,:) = T_vap_out;
    Quality(i) = 0;
    i = i+1;

end

%% Condenser 

%% Downward tray

i = 1;
for quality_out = 0:0.1:1
    i
    x_in = x_data(i,:);
    y_out = y_data(i,:);
    T_liq = T_data(i,1);
    T_vap = 89.1978;
    quality = quality_data(i);
    rg_in = rho_g_out(i);
    [x_out,y_out,T_liq_out,T_vap_out] = tray_downward(y_out,x_in,T_liq,T_vap,m_liq,m_vap);
    x_top(i,:) = x_out;
    y_top(i,:) = y_out;
    T_liq_sat(i,:) = T_liq_out;
    T_vap_sat(i,:) = T_vap_out;
    Quality(i) = 0;
    i = i+1;

end

%%

clf
figure(1)
hold on 
color = ['g' 'b' 'r'];
for k=1:3
    if k==3
        plot(quality_data,10*x_top(:,k),color = color(k))
        plot(quality_data,10*y_top(:,k),'--',color = color(k) )
    else
        plot(quality_data,x_top(:,k),color = color(k))
        plot(quality_data,y_top(:,k),'--',color = color(k) )
    end

end
hold off

%%
figure(2)
plot(quality_data(1:end),T_liq_sat,'blue')
hold on
plot(quality_data(1:end),T_vap_sat,'black')
ylabel('Temperature (K)')
xlabel('Reboiler Outlet Quality (mass)')
plotfixer


