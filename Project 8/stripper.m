function [x_in,out,tray,feed] = stripper(x_out,quality,N_trays)
% Given a guess on x_out (liquid product), number of trays and quality of
% the reboiler
% Use re boiler and go up the trays 
% Returns the liquid composition at the top
global P
P = 1e5;

% Reboiler 
[T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, y_out, rho_f_out, rho_g_out, Q_req]=reboiler_bwd(x_out, quality);
y_in = y_out;
T_liq = T_in;
T_vap = T_out;
m_vap = quality;
m_liq = 1;
x_in;
out.T = T_out;
out.x = x_out(1);
out.y = y_out(1);
out.N2 = quality*out.y+(1-quality)*out.x;
out.xS = s_crT(x_out,rl_cTP(x_out,T_out,P),T_out);
out.yS = s_crT(y_out,rv_cTP(y_out,T_out,P),T_out);
out.N2s = quality*out.yS+(1-quality)*out.xS;
out_rate = (1-quality)/M_c(x_out);
out.liq = 1;
out.vap = 0;
out.tot = out.liq + out.vap;

% Trays
for k=1:N_trays    
    tray(k).T = T_liq;
    tray(k).x = x_in(1);
    tray(k).xS = s_crT(x_in,rl_cTP(x_in,T_liq,P),T_liq);
    tray(k).vap = m_vap/M_c(y_in)/out_rate;

    [x_out,y_out,T_liq_out,T_vap_out,m_vap_out,m_liq_in] = tray_upward(x_in,y_in,T_liq,T_vap,m_vap,m_liq);
    x_in = x_out;
    y_in = y_out;
    T_liq = T_liq_out;
    T_vap = T_vap_out;
    m_liq = m_liq_in;
    m_vap = m_vap_out;
    quality = m_vap_out/(m_vap_out+m_liq_in);

    
    tray(k).y = y_in(1);
    tray(k).N2 = quality*y_in(1)+(1-quality)*tray(k).x;
    tray(k).yS = s_crT(y_in,rv_cTP(y_in,T_liq,P),T_liq);
    tray(k).N2s = quality*tray(k).yS+(1-quality)*tray(k).xS;
    tray(k).liq = m_liq_in/M_c(x_in)/out_rate;
    tray(k).tot = tray(k).vap + tray(k).liq;
end

feed.liq = m_liq_in/M_c(x_in)/out_rate;
feed.vap = m_vap_out/M_c(y_in)/out_rate;
feed.tot = feed.liq+feed.vap;