function left_spin = Liquid_Spinodal_iT(temp)
    rmax = 100;
    rmin = 0.1;
    dr = 0.001;

    cutoff = 0.001;

    % Preallocate storage...
    Prisotherm  = zeros((rmax-rmin)/dr);
    risotherm   = zeros((rmax-rmin)/dr);
    T = temp;
    i = 1;
    for r=rmin:dr:rmax+dr
        r;
        Prisotherm(i) = P_irT(ispecies,r,T);
        risotherm(i) = r;
        i = i+1;
    end
    
    left_spin = 0;

    for j = 1:length(Prisotherm)-1
        if abs(Prisotherm(j+1)-Prisotherm(j)) <=cutoff
            left_spin = risotherm(j);
            break
        end
    end

end