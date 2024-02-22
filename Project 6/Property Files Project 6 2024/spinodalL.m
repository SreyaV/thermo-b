function left_spin = spinodalL(ispecies,temp)
    rmax = 100;
    rmin = 0.1;
    dr = 0.01;

    cutoff = 300;

    % Preallocate storage...
    dPrisotherm  = zeros(1,(rmax-rmin)/dr);
    risotherm   = zeros(1,(rmax-rmin)/dr);

    T = temp;
    i = 1;
    for r=rmin:dr:rmax+dr
        dPrisotherm(1,i) = dPdr_irT(ispecies,r,T);
        risotherm(1,i) = r;
        i = i+1;
    end
    
    left_spin = 0;

    for j = 1:length(dPrisotherm)-1
        if abs(dPrisotherm(1,j)) <=cutoff
            left_spin = risotherm(j);
            break
        end
    end

end