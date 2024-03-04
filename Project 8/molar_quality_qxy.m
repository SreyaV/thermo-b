function q_mole = molar_quality_qxy(q_mass,x,y)
% Return the molar quality given the mass quality and compositions.
Mg = M_c(y);
Mf = M_c(x);
q_mole = q_mass/Mg/(q_mass/Mg+(1-q_mass)/Mf);