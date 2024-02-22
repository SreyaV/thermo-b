function right_spin = spinodalV(ispecies,temp)
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
        r;
        dPrisotherm(1,i) = dPdr_irT(ispecies,r,T);
        risotherm(1,i) = r;
        i = i+1;
    end
    
    right_spin = 0;

    for j = length(dPrisotherm)-1:-1:1
        if abs(dPrisotherm(1,j)) <=cutoff
            right_spin = risotherm(1,j);
            break
        end
    end

end