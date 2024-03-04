function q_mass = mass_quality_qxy(q_mole,x,y)
% Return the mass quality given the molar quality and compositions.
Mg = M_c(y);
Mf = M_c(x);
q_mass = (q_mole*Mg)/(q_mole*Mg+(1-q_mole)*Mf);