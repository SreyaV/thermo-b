%% T-s diagramm nH2
clear all 

Setup_Props_i;

ispecies = nH2;
Tcrit = Tcrit_i(ispecies);
Ttrip = Ttrip_i(ispecies);
Pcrit = Pcrit_i(ispecies);
Ptrip = Ptrip_i(ispecies);


%% Vapor Dome

dP = 1e5;
Plist = Ptrip:dP:0.98*Pcrit;
N = length(Plist);

for i = 1:length(Plist)
    P = Plist(i)
    [T rf rg] = Saturation_iP(ispecies,P);
    Tsat_liq(i)    = T;
    Tsat_vap(N-i+1)  = T;
    entropy_liq(i)   = s_irT(ispecies,rf,T)/1e3;
    entropy_vap(N-i+1) = s_irT(ispecies,rg,T)/1e3;
end

entropy = [entropy_liq s_irT(ispecies,rcrit_i(ispecies),Tcrit)/1e3 entropy_vap];
Tsat    = [Tsat_liq Tcrit Tsat_vap];

%% Triple line

steps = 100;
rftrip = rftrip_i(ispecies);
rgtrip = rgtrip_i(ispecies);
i = 1;

for x = 0:1/steps:1
    
    r = x*rftrip + (1-x)*rgtrip;
    s_trip(i) = s_irT(ispecies,r,Ttrip)/1e3;
    T_trip(i) = Ttrip;
    i = i+1;

end

%% Isobars
clear Tlist
clear slist

fprintf("Generating isobars (P in bars)...\n")

Plist = [500 200 100 50 20 10 5 2 1 .5 .2 .1]; %[bars]
steps = 20; %For T
count = 0;

for i=1:length(Plist)
    
    P = Plist(i)
    P = P*1e5;
    
    % Find saturation T
    if P<Pcrit
        [Tlower rf rg]= Saturation_iP(ispecies,P);
    else
        Tlower = Ttrip+0.1;
    end

    % Isobar around dome
    Tlist(i,:) = Tlower:(300-Tlower)/steps:300;

    for j=1:length(Tlist(i,:))
        
        T = Tlist(i,j);

        if P<Pcrit
            r = rv_iTP(ispecies,T,P);
        else
            r = rl_iTP(ispecies,T,P);
        end

        slist(i,j) = s_irT(ispecies,r,T)/1e3;

    end

    % Isobars inside the dome
    if P<Pcrit
        
        count = count +1;
        sliq =s_irT(ispecies,rf,Tlower)/1e3;
        svap = slist(i,1);

        for j = 1:length(Tlist(i,:))

            Tlist(length(Plist)+count,j) = Tlower;
            slist(length(Plist)+count,j) = (svap-sliq)/steps*(j-1)+sliq;
        end
    end
end

%% Isenthalps

hlist = [ 0.36 0.626 0.9 1.16 1.43 1.7 1.96 2.23 2.5 2.76 3.03 3.3 3.57 3.83 4.1 4.37]; %[MJ/kg]

fprintf("Generating isenthalps...\n")
steps = 10; 

slist_h = zeros(length(hlist),steps);
Tlist_h = zeros(length(hlist),steps);

for i=1:length(hlist)
    
    h = hlist(i)
    h = h*1e6;
    
    for j=1:length(Plist)
      
        P = Plist(j)*1e5;
        
        % Isenthalps around dome
        try  
            [r T rf rg] = rT_ihP(ispecies,h,P);
            slist_h(i,j) = s_irT(ispecies,r,T);
            Tlist_h(i,j) = T;

        % Isenthalps inside the dome
        catch
            disp('Error')
            Tlist_h(i,j) = Ttrip;
            r = rl_iTP(ispecies,Ttrip,P);
            slist_h(i,j) = s_irT(ispecies,r,Ttrip);

        end
    end
    
end

%% Plots
clf

% Triple Line
plot(s_trip,T_trip,'b')
hold on

% Isobars
for i=1:length(Plist)+count
    plot(slist(i,:),Tlist(i,:),'g')
end

%Isenthalps
for i=1:length(hlist)
    plot(slist_h(i,:)/1e3,Tlist_h(i,:),'r')
end

% Vapor Dome
plot(entropy,Tsat,'b')
hold off

xlabel("Specific entropy (kJ/kg-K)")
ylabel("Temperature (K)")
plotfixer


















