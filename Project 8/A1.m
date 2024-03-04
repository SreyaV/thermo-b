%% Setup
% Output reboiler- temperature, composition, density (liquid, vapour), complementary outlet phase, feed-specific heat transfer required 
clc;clear all;close all;
global N2 O2 Ar N P
N=3;
Setup_Air_Props;
P=100000;
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
%%
plotA1();