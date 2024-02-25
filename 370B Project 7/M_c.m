function MM = M_c(c)
% Return the molecular mass for mixture c.
% C.F. Edwards, 2-16-10

global M_i

MM = 0;
for i=1:1:length(c)
    MM = MM + c(i)*M_i(i);
end