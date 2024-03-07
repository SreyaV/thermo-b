function [x_out,y_out,T_liq_out,T_vap_out] = tray_downward(y_out,x_in,T_liq,T_vap,m_liq,m_vap)
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
    [T_liq_out rf_out rg_in x_out] = Fast_Dew_cP(y_out,P);
    
    % Above quantities are known (quality is in mass here)
    n_top_vap = y_out*quality   /M_c(y_out); %moles
    n_top_liq = x_in*(1-quality)/M_c(x_in); %moles

    % Iterate on mass flow rate of liquid and check using species continuity
    % and energy conservation
    % Method used : Bisection

    m_inf = 0;
    m_up  = 10;
    
    delta = 5e5;
    while delta > 1e3
        
        m = (m_inf + m_up)/2 % Flow rate liquid (downwards)
        N = m/M_c(x_out); % Molar flow rate 

        % Species continuity (in moles)
        n_bot_vap = n_top_vap + N*x_out - n_top_liq; %moles

        if n_bot_vap > 0
        
            % Find temperature of the vapour below
            N_bot = sum(n_bot_vap);
            y_in = n_bot_vap/sum(n_bot_vap);
            [T_vap_out rf rg_in y] = Fast_Dew_cP(y_in,P);
            T_vap_out

            % Energy conservation (in mass)
            h_bot_liq =  m                 * h_crT(x_out , rf_out, T_liq_out);
            h_bot_vap =  N_bot*M_c(y_in)   * h_crT(y_in , rg_in, T_vap_out);
            h_top_liq =  m_liq             * h_crT(x_in, rl_cTP(x_in,T_liq,P), T_liq);
            h_top_vap =  m_vap             * h_crT(y_out, rg_out, T_vap);

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