function [x_out,y_out, T_liq_out,T_vap_out, m_vap_in, m_liq_out] = tray_downward(x_in, y_in, T_liq,T_vap,m_liq,m_vap)
% Ideal Tray (downwards solving)
% Input : Liquid composition from above x_in, Vapour composition emerging
% from the tray y_out and associated Temperatures
% Output : Vapour from below y_in, Liquid composition from
% below and temperatures

global P

if T_liq>T_vap

    disp('Error in solving the tray upwards')

else 
    
    % Find Temperature and Composition of liquid phase (below)
    [T_liq_out rg_in rf_out x_out] = Fast_Dew_cP(y_in,P);
    
    % Above quantities are known (quality is in mass here)
    n_top_vap = y_in*m_vap   /M_c(y_in); %moles
    n_top_liq = x_in*m_liq/M_c(x_in); %moles

    % Iterate on mass flow rate of liquid and check using species continuity
    % and energy conservation
    % Method used : Bisection

    m_inf = 0;
    m_up  = 4;
    
    delta = 5e5;
    while delta > 1e3
        
        m = (m_inf + m_up)/2 % Flow rate liquid (downwards)
        N = m/M_c(x_out); % Molar flow rate 

        % Species continuity (in moles)
        n_bot_vap = n_top_vap + N*x_out - n_top_liq; %moles

        if n_bot_vap > 0
        
            % Find temperature of the vapour below
            N_bot = sum(n_bot_vap);
            y_out = n_bot_vap/sum(n_bot_vap);
            [T_vap_out rg_out rf_in y] = Fast_Dew_cP(y_out,P);
            T_vap_out
            mass_check = M_c(x_out)*N_bot + m_vap - m_liq -m

            % Energy conservation (in mass)
            h_bot_liq =  m                 * h_crT(x_out , rf_out, T_liq_out);
            h_bot_vap =  N_bot*M_c(y_in)   * h_crT(y_in , rg_out, T_vap_out);
            h_top_liq =  m_liq             * h_crT(x_in, rl_cTP(x_in,T_liq,P), T_liq);
            h_top_vap =  m_vap             * h_crT(y_out, rg_in, T_vap);

      
            delta = h_top_liq + h_bot_vap - h_bot_liq -h_top_vap

            if delta < 0
                m_inf = m;
            else
                m_up = m;      
            end
       
        else

            m_inf = m;
        
        end





        delta = abs(delta);

    end
    
    m_vap_in = M_c(y_out) * N_bot;
    m_liq_out = m;
    %above up is in
    %below down is out
    %liquid is out
    %vapor is in

end