function [x_in] = stripper(x_out,quality,N_trays)
% Given a guess on x_out (liquid product), number of trays and quality of
% the reboiler
% Use re boiler and go up the trays 
% Returns the liquid composition at the top

% Reboiler 
[T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, y_out, rho_f_out, rho_g_out, Q_req]=reboiler_bwd(x_out, quality);
y_in = y_out;
T_liq = T_in;
T_vap = T_out;
m_vap = quality;
m_liq = 1;
x_in;

% Trays
for k=1:N_trays    
    [x_out,y_out,T_liq_out,T_vap_out,m_vap_out,m_liq_in] = tray_upward(x_in,y_in,T_liq,T_vap,m_vap,m_liq);
    x_in = x_out;
    y_in = y_out;
    T_liq = T_liq_out;
    T_vap = T_vap_out;
    m_liq = m_liq_in;
    m_vap = m_vap_out;
end