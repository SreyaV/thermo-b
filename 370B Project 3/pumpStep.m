function [Pret,Tret,finalState,pathStates] = pump(fluid,P1,P2,eta_p,nsteps)
    vaporFraction(fluid)
    try
        set(fluid,'T',T1,'P',P1);
    catch % if the fluid is in its saturated state, set using its vapor fraction
        set(fluid,'P',P1,'Vapor',vaporFraction(fluid));
    end

    dP = (P2 - P1)/nsteps;
    Pret = P1:dP:P2;
    Tret = 0*Pret;
    hret = 0*Pret;
    sret = 0*Pret;
    vF = 0*Pret; 

    sw1 = entropy_mass(fluid);
    hw1 = enthalpy_mass(fluid);
    s   = sw1;
    h   = hw1;
    for iter = 1:1:length(Pret)
        % Find the isentropic state for the new pressure.
        set(fluid,'S',s,'P',Pret(iter));
        % Get the isentropic enthalpy difference.
        hs = enthalpy_mass(fluid);
        dhs = hs - h; 
        % The actual difference is larger due to inefficiency.
        h = h + dhs/eta_p;
        % Find the new state.
        set(fluid,'H',h,'P',Pret(iter));
        % Get the entropy for the next step in the process.
        s = entropy_mass(fluid);
        Tret(iter) = temperature(fluid);
        hret(iter) = enthalpy_mass(fluid);
        sret(iter) = entropy_mass(fluid);
        % if the fluid is a fluid, get its vapor fraction
        try
            vF(iter) = vaporFraction(fluid);
        catch
            vF(iter) = nan;
        end

    end
    
    finalState.T = Tret(end);
    finalState.P = Pret(end);
    finalState.H = enthalpy_mass(fluid);
    finalState.S = entropy_mass(fluid);
    
    pathStates.T = Tret;
    pathStates.P = Pret;
    pathStates.H = hret;
    pathStates.S = sret;
    pathStates.vF = vF;
    
    % find out exergy if required else forget about it
    pathStates.X = 0;
    pathStates.Xf = 0;

end
