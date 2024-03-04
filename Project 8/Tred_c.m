function T = Tred_c(c)
% Find the reducing T for this composition.
% C.F. Edwards, 2-16-10

global Tcrit_i
global Zeta

T = 0;
for i=1:1:length(c)
    T = T + c(i)*Tcrit_i(i);
end
for i=1:1:length(c)-1
    for j=i+1:1:length(c)
        T = T + c(i)*c(j)*Zeta(i,j);
    end
end