function [P_ideal_bubble,rho_f_ideal,rho_g_ideal,x_ideal_bubble]=i_Bubble(X,T,varargin)
% Check CFE's Bubble_cT.m for better documentation lol
global toler N2 O2 Ar Ru M_i
N=length(X);
x_ideal_bubble=zeros(N,1);
x_high=zeros(N,1);
x_low=zeros(N,1);
P_g_ideal=zeros(N,1);
for i=1:1:N
    if(Vapor_Spinodal_iT(i,T)==0)
        P_g_ideal(i)=P_crT(X,varargin{2},varargin{1});
    else
        P_g_ideal(i)=P_irT(i,Vapor_Spinodal_iT(i,T),T);
    end
end
Pmax=min(P_g_ideal)-(min(P_g_ideal/1000));
if(P_crT(X,Liquid_Spinodal_cT(X,T),T)<0)
    Pmin=1;
else
    Pmin=P_crT(X,Liquid_Spinodal_cT(X,T),T);
end
P_ideal_bubble=(Pmin+Pmax)/2;
rho_l_air=rl_cTP(X,T,P_ideal_bubble);
rho_v_N2=rv_iTP(N2,T,P_ideal_bubble);
rho_v_O2=rv_iTP(O2,T,P_ideal_bubble);
rho_v_Ar=rv_iTP(Ar,T,P_ideal_bubble);  
for i=1:1:1000 
    rho_l_air=rl_cTP(X,T,P_ideal_bubble,rho_l_air);
    rho_v_N2=rv_iTP(N2,T,P_ideal_bubble,rho_v_N2);
    x_ideal_bubble(N2)=exp((mui_icrT(N2,X,rho_l_air,T)-mu_irT(N2,rho_v_N2,T))/(Ru*T));
    rho_v_O2=rv_iTP(O2,T,P_ideal_bubble,rho_v_O2);
    x_ideal_bubble(O2)=exp((mui_icrT(O2,X,rho_l_air,T)-mu_irT(O2,rho_v_O2,T))/(Ru*T));
    rho_v_Ar=rv_iTP(Ar,T,P_ideal_bubble,rho_v_Ar);
    x_ideal_bubble(Ar)=exp((mui_icrT(Ar,X,rho_l_air,T)-mu_irT(Ar,rho_v_Ar,T))/(Ru*T));
    mass_mix=(x_ideal_bubble(N2)*M_i(N2))+(x_ideal_bubble(O2)*M_i(O2))+x_ideal_bubble(Ar)*M_i(Ar);
    vol_mix=(x_ideal_bubble(N2)*M_i(N2)/rho_v_N2)+(x_ideal_bubble(O2)*M_i(O2)/rho_v_O2)+(x_ideal_bubble(Ar)*M_i(Ar)/rho_v_Ar);
    rho_v_mix=mass_mix/vol_mix;
    if((abs(sum(x_ideal_bubble)-1)<toler))
        rho_g_ideal=rho_v_mix;
        rho_f_ideal=rho_l_air;
        x_ideal_bubble=real(x_ideal_bubble/sum(x_ideal_bubble));
        return
    end

    % 2nd order Newton-Raphson
    rho_l_air=rl_cTP(X,T,P_ideal_bubble-5,rho_l_air);
    rho_v_N2=rv_iTP(N2,T,P_ideal_bubble-5,rho_v_N2);
    x_low(N2)=exp((mui_icrT(N2,X,rho_l_air,T)-mu_irT(N2,rho_v_N2,T))/(Ru*T));
    rho_v_O2=rv_iTP(O2,T,P_ideal_bubble-5,rho_v_O2);
    x_low(O2)=exp((mui_icrT(O2,X,rho_l_air,T)-mu_irT(O2,rho_v_O2,T))/(Ru*T));
    rho_v_Ar=rv_iTP(Ar,T,P_ideal_bubble-5,rho_v_Ar);
    x_low(Ar)=exp((mui_icrT(Ar,X,rho_l_air,T)-mu_irT(Ar,rho_v_Ar,T))/(Ru*T));
    rho_l_air=rl_cTP(X,T,P_ideal_bubble+5,rho_l_air);
    rho_v_N2=rv_iTP(N2,T,P_ideal_bubble+5,rho_v_N2);
    x_high(N2)=exp((mui_icrT(N2,X,rho_l_air,T)-mu_irT(N2,rho_v_N2,T))/(Ru*T));
    rho_v_O2=rv_iTP(O2,T,P_ideal_bubble+5,rho_v_O2);
    x_high(O2)=exp((mui_icrT(O2,X,rho_l_air,T)-mu_irT(O2,rho_v_O2,T))/(Ru*T));
    rho_v_Ar=rv_iTP(Ar,T,P_ideal_bubble+5,rho_v_Ar);
    x_high(Ar)=exp((mui_icrT(Ar,X,rho_l_air,T)-mu_irT(Ar,rho_v_Ar,T))/(Ru*T));


    delF=(sum(x_high)-sum(x_low))/(2*5);                                % derivative
    P_ideal_bubble=P_ideal_bubble-(sum(x_ideal_bubble)-1)/delF;         % New Pressure estimate
    if(P_ideal_bubble<Pmin)
        P_ideal_bubble=Pmin;
    end
    if(Pmax<P_ideal_bubble)
        P_ideal_bubble=Pmax;
    end
end
P_ideal_bubble=0;                       % fallback to 0s
rho_f_ideal=0;
rho_g_ideal=0;
x_ideal_bubble=zeros(1,N);