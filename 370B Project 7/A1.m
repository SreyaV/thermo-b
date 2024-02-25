clc;clear;close;
%% Setup EOS Parameters from https://srd.nist.gov/jpcrdreprint/1.1285884.pdf
M_air=0.02896546*1000;
T_j=132.6312;
P_j=3.78502e6;
rho_j=10.4477*M_air;
T_p=132.6035;
P_p=3.7891e6;
rho_p=11.0948*M_air;
T_c=132.5306;
P_c=3.7860e6;
rho_c=11.8308*M_air;
N=3;
Setup_Air_Props;
X=zeros(N,1);
X([N2, O2, Ar]) = [0.7812, 0.2096, 0.0092];
[t,r]=Pr_Inflection_c(X);
p1=P_crT(X,r,t);
%% Dew

for T=60:1:t
    [P_ideal,rho_g_ideal,rho_f_ideal,x_ideal]=i_Dew(X,T,t,r);
    if(P_ideal==0)
        break
    end
    dew_data(T-59, :)=[T,P_ideal,rho_g_ideal,rho_f_ideal];
    x_ideal_dew(T-59,:)=x_ideal;
    [P_real,rho_g_real,rho_f_real,x_real]=Dew_cT(X,T,P_ideal,x_ideal,rho_f_ideal,rho_g_ideal);
    dew_data_real(T-59, :)=[T,P_real,rho_g_real,rho_f_real];
    x_real_dew(T-59, :)=x_real;
end

%% Bubble

for T=60:1:t
    [P_ideal_bubble,rho_f_ideal_bubble,rho_g_ideal_bubble,x_ideal_bubble]=i_Bubble(X,T,t,r);
    if(P_ideal_bubble==0)
        break
    end
    bubble_data(T-59, :)=[T, P_ideal_bubble, rho_f_ideal_bubble, rho_g_ideal_bubble];
    x_ibubble(T-59,:)=x_ideal_bubble;
    [P_real,rho_f_real,rho_g_real,x_real]=Bubble_cT(X,T,P_ideal_bubble,x_ideal_bubble,rho_f_ideal_bubble,rho_g_ideal_bubble);
    bubble_data_real(T-59,:)=[T,P_real,rho_f_real,rho_g_real];
    x_ibubble_real(T-59,:) = x_real;
end
plotA1
