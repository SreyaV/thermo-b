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
    x_data(i, :)=x_in;
    y_data(i, :)=y_out;
    T_data(i, 1:2)=[T_in, T_out];
    Q_data(i)=Q_req;
    quality_data(i)=quality_out;
    i=i+1;
end

%% Upward tray

i = 9;
for quality_out = 0.8%:0.1:1
    i
    x_in = x_data(i,:);
    y_in = y_data(i,:);
    T_vap = 89.5;
    T_liq = T_data(i,1);
    quality = quality_data(i);
    [x_out,y_out,T_liq_out,T_vap_out] = tray_upward(x_in,y_in,T_liq,T_vap,quality);
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
plot(quality_data(1:end-3),T_liq_sat)
hold on
plot(quality_data(1:end-3),T_vap_sat)
plotfixer


