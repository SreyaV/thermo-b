function [x_out,y_out,T_liq_out,T_vap_out,m_vap_out,m_liq_in] = tray_upward(x_in,y_in,T_liq,T_vap,m_vap,m_liq)
% Ideal Tray (upwards solving)
% Input : Vapour composition from below y_in, Liquid composition emerging
% from the tray x_in and associated Temperatures
% Output : Vapour emerging from the tray y_out, Liquid composition from
% above and temperatures

global P

if T_liq>T_vap

    disp('Error in solving the tray upwards')

else 
    
    % Find Temperature and Composition of vapour phase (for an ideal tray)
    [T_vap_out rf_in rg_out y_out] = Fast_Bubble_cP(x_in,P);
    
    % Bottom quantities are known (quality is in mass here)
    n_bot_vap = y_in*m_vap /M_c(y_in); %moles
    n_bot_liq = x_in*m_liq /M_c(x_in); %moles

    % Iterate on mass flow rate of vapor and check using species continuity
    % and energy conservation
    % Method used : Bisection

    m_inf = 0;
    m_up  = 10;
    
    delta = 5e5;
    while delta > 1e2
        
        m = (m_inf + m_up)/2; % Flow rate vapor (upwards)
        N = m/M_c(y_out); % Molar flow rate 

        % Species continuity (in moles)
        n_top_liq = n_bot_liq + N*y_out -n_bot_vap; %moles

        if n_top_liq > 0
        
            % Find temperature of the liquid above
            N_top = sum(n_top_liq);
            x_out = n_top_liq/N_top;
            [T_liq_out rf_out rg y] = Fast_Bubble_cP(x_out,P);
            mass_check = M_c(x_out)*N_top + m_vap - m_liq -m;

            % Energy conservation (in mass)
            h_bot_liq =  m_liq              * h_crT(x_in , rl_cTP(x_in,T_liq,P), T_liq);
            h_bot_vap =  m_vap              * h_crT(y_in , rv_cTP(y_in,T_vap,P), T_vap);
            h_top_liq =  M_c(x_out) * N_top * h_crT(x_out, rl_cTP(x_out,T_liq_out,P), T_liq_out);
            h_top_vap =  m                  * h_crT(y_out, rv_cTP(y_out,T_vap_out,P), T_vap_out);

            delta = h_top_liq + h_bot_vap - h_bot_liq -h_top_vap;

            if delta > 0
                m_inf = m;
            else
                m_up = m;      
            end
       
        else

            m_inf = m;
        
        end





        delta = abs(delta);

    end
    
    m_vap_out = m;
    m_liq_in = M_c(x_out) * N_top;

end