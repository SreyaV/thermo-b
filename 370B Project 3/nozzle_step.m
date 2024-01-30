function fluid = nozzle(fluid, ratio_n)
    P = pressure(fluid);
    h = enthalpy_mass(fluid);
    set(fluid,'P',P*ratio_n,'H',h);

end
