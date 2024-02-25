function R = R_c(c)
% Return the specific gas constant for mixture c.
% C.F. Edwards, 2-16-10

global M_i Ru

MM = 0;
for i=1:1:length(c)
    MM = MM + c(i)*M_i(i);
end

R = Ru/MM;
