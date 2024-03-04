function [T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, y_out, rho_f_out, rho_g_out, Q_req]=reboiler_bwd(x_out, quality)
% Output- temperature, composition, density (liquid, vapour), complementary outlet phase, feed-specific heat transfer required 
global N2 O2 Ar P
[T_out, rho_f_out, rho_g_out, y_out]=Bubble_cP(x_out, P);
M=M_c(y_out)*M_c(x_out)/(quality*M_c(x_out)+M_c(y_out)-quality*M_c(y_out));
x_in(N2)=y_out(N2)*quality*M/M_c(y_out)+x_out(N2)-quality*M*x_out(N2)/M_c(y_out);
x_in(O2)=y_out(O2)*quality*M/M_c(y_out)+x_out(O2)-quality*M*x_out(O2)/M_c(y_out);
x_in(Ar)=y_out(Ar)*quality*M/M_c(y_out)+x_out(Ar)-quality*M*x_out(Ar)/M_c(y_out);
[T_in, rho_f_in, rho_g_in, y_in]=Bubble_cP(x_in, P);
Q_req=quality*h_crT(y_out,rho_g_out,T_out)-quality*h_crT(x_out,rho_f_out,T_out)+h_crT(x_out,rho_f_out,T_out)-h_crT(x_in,rho_f_in,T_in);
