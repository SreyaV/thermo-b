clc; clear; close all;

%% Setup Anode Gas Humidified Hydrogen
gas=Solution('GRI30.yaml');
iso_temp=1273.15;
GDL_pressure=100000;
X_A=zeros(nSpecies(gas), 1);
X_A(speciesIndex(gas, 'H2O'))=0.03;          
X_A(speciesIndex(gas, 'H2'))=0.97;
set(gas, 'T', iso_temp, 'P', GDL_pressure, 'X', X_A);
mu_A=chemPotentials(gas)/1000;

%% Setup Cathode Gas Engineering Air
X_C=zeros(nSpecies(gas), 1);
X_C(speciesIndex(gas, 'N2'))=0.79;
X_C(speciesIndex(gas, 'O2'))=0.21;
set(gas,'T', iso_temp,'P', GDL_pressure,'X', X_C);
mu_C=chemPotentials(gas)/1000;
Affinity0=(mu_A(speciesIndex(gas, 'H2'))-mu_A(speciesIndex(gas, 'H2O'))+0.5*mu_C(speciesIndex(gas, 'O2')));

%% With Loss at reference T=1000C
R=8.31447;
F=96.5*10^3;
a=3.640*10^-4;
b=2.334;
m_H2=2;
m_H2O=18;
P_H2=12.9583;
P_H2O=217.6659;
T_H2=33.19;
T_H2O=647.1;
diffusion_H2_H2O=(((P_H2*P_H2O)^(0.33)*(T_H2*T_H2O)^(0.4167)*(1/m_H2+1/m_H2O)^(0.5)*a*(iso_temp/sqrt(T_H2*T_H2O))^b)/(GDL_pressure/101325))/(10000);
a=2.745*10^-4;
b=1.823;
m_O2=32;
m_N2=28;
P_O2=49.77;
P_N2=33.55;
T_O2=154.6;
T_N2=126.2;
diffusion_O2_N2=(((P_O2*P_N2)^(0.33)*(T_O2*T_N2)^(0.4167)*(1/m_O2+1/m_N2)^(0.5)*a*(iso_temp/sqrt(T_O2*T_N2))^b)/(GDL_pressure/101325))/(10000);

%% Temp Sweep b/w 600-1300C
k0=15;
Activation_Energy_A=F;
exchange_density_C=1000;
T0=1273.15;
Activation_A=1000*100;
i=1;
for Temp=600:10:1300
    T_data(i)=Temp;
    T=Temp+273.15;
    k_T(i)=k0*exp(-(Activation_Energy_A/R)*((1/T)-(1/T0)));
    a=3.640*10^-4;
    b=2.334;
    m_H2=2;
    m_H2O=18;
    P_H2=12.9583;
    P_H2O=217.6659;
    T_H2=33.19;
    T_H2O=647.1;
    diffusion_data(1, i)=(((P_H2*P_H2O)^(0.33)*(T_H2*T_H2O)^(0.4167)*(1/m_H2+1/m_H2O)^(0.5)*a*(T/sqrt(T_H2*T_H2O))^b)/(GDL_pressure/101325))/(10000); 
    a=2.745*10^-4;
    b=1.823;
    m_O2=32;
    m_N2=28;
    P_O2=49.77;
    P_N2=33.55;
    T_O2=154.6;
    T_N2=126.2;
    diffusion_data(2, i)=(((P_O2*P_N2)^(0.33)*(T_O2*T_N2)^(0.4167)*(1/m_O2+1/m_N2)^(0.5)*a*(T/sqrt(T_O2*T_N2))^b)/(GDL_pressure/101325))/(10000);
    set(gas, 'T', T, 'P', GDL_pressure, 'X', X_A);
    mu_A=chemPotentials(gas)/1000;
    set(gas, 'T', T, 'P', GDL_pressure, 'X', X_C);
    mu_C=chemPotentials(gas)/1000;
    Affinity_data(i)=mu_A(speciesIndex(gas, 'H2'))-mu_A(speciesIndex(gas, 'H2O'))+0.5*mu_C(speciesIndex(gas, 'O2'));
    current_density_C_data(i)=exchange_density_C*exp(-(Activation_A/R)*((1/T)-(1/T0)));
    i=i+1;
end

%% Plot
hold on
plot(T_data,Affinity_data/Affinity0)
plot(T_data,diffusion_data(2, :)/diffusion_O2_N2)
plot(T_data,diffusion_data(1, :)/diffusion_H2_H2O)
plot(T_data,k_T/k0)
plot(T_data,current_density_C_data/exchange_density_C)
hold off
ylabel('Normalized data (T_{ref}=1000\circC)')
xlabel('T\circC')
legend('Overall Reaction Affinity', 'O_2/N_2 Diffusivity', 'H_2/H_2O Diffusivity', 'Ionic Conductivity of YSZ', 'Cathode Exchange Current Density')
plotfixer
