function right_spin = Vapor_Spinodal_iT(temp)
    rmax = 100;
    rmin = 0.1;
    dr = 0.001;

    cutoff = 0.0001;

    % Preallocate storage...
    dPrisotherm  = zeros((rmax-rmin)/dr);
    risotherm   = zeros((rmax-rmin)/dr);
    T = temp;
    i = 1;
    for r=rmin:dr:rmax+dr
        r;
        dPrisotherm(i) = dPdr_irT(ispecies,r,T);
        risotherm(i) = r;
        i = i+1;
    end
    
    right_spin = 0;

    for j = length(dPrisotherm)-1:-1:1
        if abs(dPrisotherm(j)) <=cutoff
            right_spin = risotherm(j);
            break
        end
    end

end