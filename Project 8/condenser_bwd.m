function [T_in, y_in, x_in, rho_f_in, rho_g_in, T_out, x_out, rho_f_out, rho_g_out, Q_req]=condenser_bwd(y_out, quality)
% Output- temperature, composition, density (liquid, vapour), complementary outlet phase, out-specific heat transfer required 
global N2 O2 Ar P
[T_out, rho_g_out, rho_f_out, x_out]=Fast_Dew_cP(y_out, P); % Dew point calculations for the out
M=M_c(x_out)*M_c(y_out)/(quality*M_c(y_out)+M_c(x_out)-quality*M_c(x_out));
y_in(N2)=x_out(N2)*quality*M/M_c(x_out)+y_out(N2)-quality*M*y_out(N2)/M_c(x_out);
y_in(O2)=x_out(O2)*quality*M/M_c(x_out)+y_out(O2)-quality*M*y_out(O2)/M_c(x_out);
y_in(Ar)=x_out(Ar)*quality*M/M_c(x_out)+y_out(Ar)-quality*M*y_out(Ar)/M_c(x_out);
[T_in, rho_g_in, rho_f_in, x_in]=Fast_Dew_cP(y_in, P); % inlet conditions based on the output composition
Q_req=h_crT(y_in,rho_g_in,T_in) - ( quality*h_crT(y_out,rho_g_out,T_out) + (1 - quality) * h_crT(x_out,rho_f_out,T_out) );