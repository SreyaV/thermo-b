function [tray,out,] = stripper_wrap(Reboiler_quality,N_trays)


% Set Binary air
z_feed([N2, O2])=[0.71 0.29];
z = z_feed;
[q V y_feed x_feed rg rf] = Flash_zTP(z,79.7,P);
%[T x_feed y_feed rf rg] = Flash_zqP(z,0.5,P);


% Make a guess for x_out
x_out = [0.1 0.9];

% Newton-Raphson
tol = 1e-4;
delta = 5;
iterations = 0;

while delta>tol
    
    iterations = iterations+1
    disp("Guess x_out")
    x_out

    [x_in,out,tray] = stripper(x_out,reboiler_quality,N_trays)

    % Derivative to find the new NR iterate (Euler forward)
    dx      = 0.01;
    x_in_dx = stripper(x_out+[dx -dx],reboiler_quality,N_trays);
    diff    = (x_in_dx(1)-x_in(1))/dx;
    
    % New iterate
    x_temp(1) = x_out(1) - (x_in(1)-x_feed(1))/diff;
    x_temp(2) = 1-x_out(1);

    % Keep it in range 0-1 (bisection of last interval)
    if x_temp(1)>1
        x_temp(1) = (1 + x_out(1))/2;
    end
    if x_temp(1)<0
        x_temp(1) = x_out(1)/2;
    end
    x_temp(2) = 1-x_temp(1);

    x_out(1) = x_temp(1);
    x_out(2) = x_temp(2);

    delta = norm(x_in-x_feed)

end


end


