function [T_in, x_in, y_in, rho_f_in, rho_g_in, T_out, x_out, rho_f_out, rho_g_out, Q_req]=condenser_bwd(y_out,quality)
global N2 O2 Ar P
[T_out, rho_g_out, rho_f_out, x_out]=Fast_Dew_cP(y_out, P);
M=M_c(y_out)*M_c(x_out)/(quality*M_c(x_out)+M_c(y_out)-quality*M_c(y_out));
y_in(N2)=y_out(N2)*quality*M/M_c(y_out) + x_out(N2)*(1-quality*M/M_c(y_out)); 
y_in(O2)=y_out(O2)*quality*M/M_c(y_out) + x_out(O2)*(1-quality*M/M_c(y_out)); 
y_in(Ar)=y_out(Ar)*quality*M/M_c(y_out) + x_out(Ar)*(1-quality*M/M_c(y_out)); 
[T_in, rho_g_in, rho_f_in, x_in]=Fast_Dew_cP(y_in,P);
Q_req=h_crT(y_in,rho_g_in,T_in)-h_crT(x_out,rho_f_out,T_out)-quality*h_crT(y_out,rho_g_out,T_out)+quality*h_crT(x_out,rho_f_out,T_out);
