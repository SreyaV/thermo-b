function [x_out,y_out,T_liq_out,T_vap_out] = tray_upward(x_in,y_in,T_liq,T_vap,quality)
% Ideal Tray (upwards solving)
% Input : Vapour composition from below x_in, Liquid composition emerging
% from the tray y_in and associated Temperatures
% Output : Vapour emerging from the tray x_out, Liquid composition from
% above and temperatures

global P

if T_liq>T_vap

    disp('Error in solving the tray upwards')

else 
    
    % Find Temperature and Composition of vapour phase (for an ideal tray)
    [T_vap_out rf rg y_out] = Fast_Bubble_cP(x_in,P);
    
    % Bottom quantities are known (quality is in mass here)
    n_bot_vap = y_in*quality    /M_c(y_in); %moles
    n_bot_liq = x_in*(1-quality)/M_c(x_in); %moles

    % Iterate on mass flow rate of vapor and check using species continuity
    % and energy conservation
    % Method used : Bisection

    m_inf = 0;
    m_up  = 10;
    
    delta = 5e5;
    while delta > 1e4
        
        m = (m_inf + m_up)/2; % Flow rate vapor
        N = m/M_c(y_out) % Molar flow rate 

        % Species continuity (in moles)
        n_top_liq = n_bot_liq + N*y_out -n_bot_vap; %moles

        if n_top_liq > 0
        
            % Find temperature of the liquid above
            N_top = sum(n_top_liq);
            x_out = n_top_liq/N_top
            [T_liq_out rf rg y] = Fast_Bubble_cP(x_out,P);
            T_liq_out

            % Energy conservation (in mass)
            h_bot_liq =  (1-quality)          * h_crT(x_in , rl_cTP(x_in,T_liq,P), T_liq);
            h_bot_vap =  quality              * h_crT(y_in , rv_cTP(y_in,T_vap,P), T_vap);
            h_top_liq =  M_c(x_out) * N_top   * h_crT(x_out, rl_cTP(x_out,T_liq_out,P), T_liq_out);
            h_top_vap =  M_c(y_out) * N       * h_crT(y_out, rv_cTP(y_out,T_vap_out,P), T_vap_out);

            delta = h_top_liq + h_bot_vap - h_bot_liq -h_top_vap

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


end