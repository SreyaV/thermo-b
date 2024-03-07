function [T_in, y_in, x_in, rho_f_in, rho_g_in, T_out, x_out, rho_f_out, rho_g_out, Q_req]=condenser_bwd(y_feed, quality)
% If the givens are instead x_feed and y_feed, quality can be calculated as
% the vapor fraction and passed to the condenser instead
% Output- temperature, composition, density (liquid, vapour), complementary outlet phase, feed-specific heat transfer required 
global N2 O2 Ar P
[T_in, rho_g_in, rho_f_in, x_in]=Dew_cP(y_feed, P); % Dew point calculations for the feed
M=M_c(x_in)*M_c(y_feed)/(quality*M_c(y_feed)+M_c(x_in)-quality*M_c(x_in));
y_out(N2)=x_in(N2)*quality*M/M_c(x_in)+y_feed(N2)-quality*M*y_feed(N2)/M_c(x_in);
y_out(O2)=x_in(O2)*quality*M/M_c(x_in)+y_feed(O2)-quality*M*y_feed(O2)/M_c(x_in);
y_out(Ar)=x_in(Ar)*quality*M/M_c(x_in)+y_feed(Ar)-quality*M*y_feed(Ar)/M_c(x_in);
[T_out, rho_g_out, rho_f_out, x_out]=Dew_cP(y_out, P); % outlet conditions based on the output composition
Q_req=quality*h_crT(x_in,rho_f_in,T_in)-quality*h_crT(y_feed,rho_g_in,T_in)+h_crT(y_feed,rho_g_in,T_in)-h_crT(y_out,rho_g_out,T_out);
