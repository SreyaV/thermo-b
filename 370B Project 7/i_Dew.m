function [P_ideal,rho_g_ideal,rho_f_ideal,x_ideal]=i_Dew(X,T,varargin)
% Check CFE's Dew.m for better documentation lol
global toler N2 O2 Ar Ru M_i
N=length(X);
x_ideal=zeros(N,1);
x_low=zeros(N,1);           % Below P_ideal
x_high=zeros(N,1);          % Above p_ideal
for i=1:1:N
    p_ideal(i)=P_irT(i,Liquid_Spinodal_iT(i,T),T);
end
p_ideal(N+1)=1;
Pmin=max(p_ideal)+(max(p_ideal)/1000);
Pmax=min([P_crT(X,varargin{2},varargin{1}), P_crT(X,Vapor_Spinodal_cT(X,T),T)]);
P_ideal=(Pmin+Pmax)/2;                                      % Start test in the middle


rho_air_v=rv_cTP(X,T,P_ideal);                              % initialize
rho_N2_l=rl_iTP(N2,T,P_ideal);
rho_O2_l=rl_iTP(O2,T,P_ideal);
rho_Ar_l=rl_iTP(Ar,T,P_ideal);


for i=1:1:1000                                              % Taking delta P=5 (can be changed to higher/lower
    rho_air_v=rv_cTP(X,T,P_ideal,rho_air_v);
    rho_O2_l=rl_iTP(O2,T,P_ideal,rho_O2_l);
    x_ideal(O2)=exp((mui_icrT(O2,X,rho_air_v,T)-mu_irT(O2,rho_O2_l,T))/(Ru*T));             % Difference in chem potentials is RTln(x_i)
    rho_N2_l=rl_iTP(N2,T,P_ideal,rho_N2_l);
    x_ideal(N2)=exp((mui_icrT(N2,X,rho_air_v,T)-mu_irT(N2,rho_N2_l,T))/(Ru*T));
    rho_Ar_l=rl_iTP(Ar,T,P_ideal,rho_Ar_l);
    x_ideal(Ar)=exp((mui_icrT(Ar,X,rho_air_v,T)-mu_irT(Ar,rho_Ar_l,T))/(Ru*T));
    mass_mix=(x_ideal(N2)*M_i(N2))+(x_ideal(O2)*M_i(O2))+(x_ideal(Ar)*M_i(Ar));
    v_mix=(x_ideal(N2)*M_i(N2)/rho_N2_l)+(x_ideal(O2)*M_i(O2)/rho_O2_l)+(x_ideal(Ar)*M_i(Ar)/rho_Ar_l);
    rho_l_mix=mass_mix/v_mix;                            % amagat for density


    if((abs(sum(x_ideal)-1)<toler))                     % Physically viable solution
        x_ideal=real(x_ideal/sum(x_ideal));
        rho_f_ideal=rho_l_mix;
        rho_g_ideal=rho_air_v;
        return                                          % Close the pressure iteration
    end

    % 2nd order Newton-Raphson
    rho_air_v=rv_cTP(X,T,P_ideal-5,rho_air_v);
    rho_l=rl_iTP(N2,T,P_ideal-5,rho_N2_l);
    x_low(N2)=exp((mui_icrT(N2,X,rho_air_v,T)-mu_irT(N2,rho_l,T))/(Ru*T));
    rho_l=rl_iTP(O2,T,P_ideal-5,rho_O2_l);
    x_low(O2)=exp((mui_icrT(O2,X,rho_air_v,T)-mu_irT(O2,rho_l,T))/(Ru*T));
    rho_l=rl_iTP(Ar,T,P_ideal-5,rho_Ar_l);
    x_low(Ar)=exp((mui_icrT(Ar,X,rho_air_v,T)-mu_irT(Ar,rho_l,T))/(Ru*T));


    rho_air_v=rv_cTP(X,T,P_ideal+5,rho_air_v);
    rho_l=rl_iTP(N2,T,P_ideal+5,rho_N2_l);
    x_high(N2)=exp((mui_icrT(N2,X,rho_air_v,T)-mu_irT(N2,rho_l,T))/(Ru*T));
    rho_l=rl_iTP(O2,T,P_ideal+5,rho_O2_l);
    x_high(O2)=exp((mui_icrT(O2,X,rho_air_v,T)-mu_irT(O2,rho_l,T))/(Ru*T));
    rho_l=rl_iTP(Ar,T,P_ideal+5,rho_Ar_l);
    x_high(Ar)=exp((mui_icrT(Ar,X,rho_air_v,T)-mu_irT(Ar,rho_l,T))/(Ru*T));


    delF=((sum(x_high))-(sum(x_low)))/(5);              % derivative
    P_ideal=P_ideal-((sum(x_ideal)-1)/delF);
    if(P_ideal<Pmin)
        P_ideal=Pmin;
    end
    if(Pmax<P_ideal)
        P_ideal=Pmax;
    end
end


P_ideal=0;                                          % fallback to 0s
rho_g_ideal=0;
rho_f_ideal=0;
x_ideal=zeros(1,N);