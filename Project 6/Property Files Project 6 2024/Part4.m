%% T-s diagramm nH2
clear all 

Setup_Props_i;

ispecies = nH2;
Tcrit = Tcrit_i(ispecies);
Ttrip = Ttrip_i(ispecies);
Pcrit = Pcrit_i(ispecies);
Ptrip = Ptrip_i(ispecies);


%% Vapor Dome

disp('Generating Vapor Dome for P (Pa)...\n')
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
    s_trip(i) = max(s_irT(ispecies,r,Ttrip)/1e3,s_irT(ispecies,rftrip,Ttrip)/1e3);
    T_trip(i) = Ttrip;
    i = i+1;

end

%% Isobars
clear Tlist
clear slist

fprintf("Generating isobars (P in bars)...\n")

Plist = [500 200 100 50 20 10 5 2 1 .5 .2 .1]; %[bars]
steps = 100; %For T
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
Phlist = [500 350 200 150 100 75 50 35 20 15 10 7.5 5 3.5 2 1.5 1 0.75 .5 0.35 .2 0.15 .1]; %[bars]
steps = length(Phlist);

fprintf("Generating isenthalps...(MJ/kg) \n")

slist_h = zeros(length(hlist),steps);
Tlist_h = zeros(length(hlist),steps);

for i=1:length(hlist)
    
    h = hlist(i)
    h = h*1e6;
    
    for j=1:steps
      
        P = Phlist(j)*1e5
        
        try  

            %Isenthalps around the dome
            [r T rf rg] = rT_ihP(ispecies,h,P);
            slist_h(i,j) = s_irT(ispecies,r,T)/1e3;
            Tlist_h(i,j) = T;

            %Isenthalps inside the dome
            hcrit = 0.8*1e6;
            if P<Pcrit_i(ispecies) && h<hcrit
                
                disp('Check dome')
                % Find Tsat(P) as well as sf and sg
                [Tsat_h rfsat rgsat] = Saturation_iP(ispecies,P);
                sf = s_irT(ispecies,rfsat,Tsat_h)/1e3;
                sg = s_irT(ispecies,rgsat,Tsat_h)/1e3;

                % Check if it is inside the dome
                if (sf <= slist_h(i,j)) && (slist_h(i,j)<=sg)
                    
                    disp('inside the dome, P(bars)=')
                    P*1e-5
                    hf = h_irT(ispecies,rfsat,Tsat_h);
                    hg = h_irT(ispecies,rgsat,Tsat_h);
                    quality = (h-hf)/(hg-hf);
                    slist_h(i,j) = sg*quality+sf*(1-quality);
                    Tlist_h(i,j) = Tsat_h;

                end
            end


        % Isenthalps Liquid side (hits triple T)
        catch
            disp('Liquid side : T = Trip')
            Tlist_h(i,j) = Ttrip;
            r = rl_iTP(ispecies,Ttrip,P);
            slist_h(i,j) = s_irT(ispecies,r,Ttrip)/1e3;

%             min = abs(h_irT(6,rl_iTP(6,Ttrip,P),Ttrip)-h) ;
%             Tlist_h(i,j) = Ttrip;
%             slist_h(i,j) = s_irT(ispecies,rl_iTP(ispecies,Ttrip,P),Ttrip)/1e3;
% 
%             for T=Ttrip:0.1:40
%                 htemp = h_irT(6,rl_iTP(6,T,P),T)*1e-6;
%                 
%                 if abs(htemp-h)<min
%                     Tlist_h(i,j) = T;
%                     r = rl_iTP(ispecies,T,P);
%                     slist_h(i,j) = s_irT(ispecies,r,T)/1e3;
%                 end
%             end


        end
    end
    
end

%% Plots
clf

% Triple Line
plot(s_trip,T_trip,'r')
hold on

% Isobars
for i=1:length(Plist)+count
    plot(slist(i,:),Tlist(i,:),'g')
end

% Add pressure labels.
steps = 100;
text(slist(1,steps)-10,Tlist(1,steps)+25,"P (bars) = ")
c=1;
for i=1:1:length(Plist)
    text(slist(i,steps),Tlist(i,steps)+25+5*sign(c),...
        num2str(Plist(i)),color='g')
    c=-c;
end

%Isenthalps
for i=1:length(hlist)
    plot(slist_h(i,:),Tlist_h(i,:),'b')
end

% Enthalpy labels
steps = length(Phlist);
text(slist_h(1,steps)-15,Tlist_h(1,steps)-12,"h (MJ/kg) = ")
for i=2:length(hlist)-1
    text(slist_h(i,steps)+1,Tlist_h(i,steps),...
        num2str(hlist(i)),color='b')
end
text(slist_h(1,steps),Tlist_h(1,steps)-12,...
        num2str(hlist(1)),color='b')

% Vapor Dome
plot(entropy,Tsat,'r')
plot(s_trip,T_trip,'r') % Triple line
hold off

xlabel("Specific entropy (kJ/kg-K)")
ylabel("Temperature (K)")
plotfixer
legend("off")



















