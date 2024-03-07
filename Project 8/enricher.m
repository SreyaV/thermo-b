function [y_in,out,tray] = enricher(y_out,quality,N_trays)
%function [condenser_out, below_up, heat_mass,trays] = enricher(y_out, quality ,N_trays)
%quality wrt mass
%condenser output is struct
%heat_mass = Q_mass/quality;
% condenser = phase;
% condenser.T = Tout;
% xout
% yout
% rfout
% rgout
% qmass
% ...

%Below up contains the state of the fluid that must be passed to the bottom
%tray. Check this state against the values of the fluid

% Condenser
[T_in, y_in, x_in, rho_f_in, rho_g_in, T_out, x_out, rho_f_out, rho_g_out, Q_req]=condenser_bwd(y_out, quality)


%x_in = x_out;
T_liq = T_in;
T_vap = T_out;
m_vap = quality;
m_liq = 1;
out.T = T_out;
out.x = x_out(1);
out.y = y_out(1);
out.N2 = quality*out.y+(1-quality)*out.x;
%out.s = 

% Trays
for k=1:N_trays    
    tray(k).T = T_vap;
    tray(k).y = y_in(1);

    [x_out,y_out,T_liq_out,T_vap_out, m_vap_in, m_liq_out] = tray_downward(y_out,x_in,T_liq,T_vap,m_liq,m_vap);

    x_in = x_out;
    y_in = y_out;
    T_liq = T_liq_out;
    T_vap = T_vap_out;
    m_liq = m_liq_out;
    m_vap = m_vap_in;
    quality = m_vap_in/(m_vap_in+m_liq_out);

    
    tray(k).x = x_in(1);
    tray(k).N2 = quality*y_in(1)+(1-quality)*tray(k).x;
end