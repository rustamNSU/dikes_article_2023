function Vg = Vg_x(dP, rhog, alpha)
    g = 9.81;
    k0 = 5.0e-12;
    j = 3.5;
    mu = 1.5e-5;

    kg = k0 * alpha.^j;
    Vg = -kg ./ mu .* (dP - rhog * g);
end