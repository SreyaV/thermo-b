function r = rred_c(c)
% Find the reducing density for this composition.
% C.F. Edwards, 2-16-10

global rcrit_i M_i
global Xi

% Work in molar volume since Xi is in those units.
vmred = 0;
for i=1:1:length(c)
    if(c(i) ~= 0)
        vmred = vmred + c(i)/(rcrit_i(i)/M_i(i));
    end
end
for i=1:1:length(c)-1
    for j=i+1:1:length(c)
        vmred = vmred + c(i)*c(j)*Xi(i,j);
    end
end

r = 1/(vmred/M_c(c));

%     if(c(i) ~= 0)
%         vred = vred + c(i)/rcrit(i);
%     end
% end
% r = 1/vred;