function Ql = get_host_rock_first_element_heating(newFrac, oldFrac, reservoir)
    fracElements = newFrac.get_fracture_elements();
    Tnew = newFrac.rockTemperature.T(fracElements, 1);
    Told = oldFrac.rockTemperature.T(fracElements, 1);
    dy = newFrac.rockTemperature.dy(1);
    rho = newFrac.rockTemperature.rho;
    C = newFrac.rockTemperature.C;

    Ql = rho * C * (Tnew - Told);
end
