function fluid = condenser(fluid, P_cond, ratio_c, TurbineExit_Q)

    set(fluid,'P',P1,'Vapor',TurbineExit_Q);

    P2 = P_cond*ratio_c;
    Q_final = 0;

    set(fluid,'P',P2,'Vapor',Q);

end
