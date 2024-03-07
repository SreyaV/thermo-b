function [feed,tray,out,x_out] = stripper_wrap(reboiler_quality,N_trays)
global N2 O2 P N
P = 1e5;

% Set Binary air
z_feed=[0.71 0.29];
z = z_feed;
T_feed = 79.7;
[q V y_feed x_feed rg rf] = Flash_zTP(z,T_feed,P);
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

    [x_in,out,tray,feed] = stripper(x_out,reboiler_quality,N_trays)

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

feed.T = T_feed;
feed.x = x_feed(1);
feed.y = y_feed(1);
feed.N2 = 0.5*feed.y+(1-0.5)*feed.x;
feed.xS = s_crT(x_feed,rl_cTP(x_feed,T_feed,P),T_feed);
feed.yS = s_crT(y_feed,rv_cTP(y_feed,T_feed,P),T_feed);
feed.N2s = 0.5*feed.yS+(1-0.5)*feed.xS;
feed.vap = feed.vap ;


end


