function solve_heat_eq_host_rock(newFrac, oldFrac, dt, reservoir, covering)
    fracElements = newFrac.get_fracture_elements();
    TInf = newFrac.rockTemperature.TInf;
    dy = newFrac.rockTemperature.dy;
    yc = newFrac.rockTemperature.yc;
    ny = newFrac.rockTemperature.ny;
    rho = newFrac.rockTemperature.rho;
    C = newFrac.rockTemperature.C;
    k = newFrac.rockTemperature.k;

    n = 3*ny - 4;
    I = zeros(1, n);
    J = zeros(1, n);
    V = zeros(1, n);
    b = zeros(1, ny);
    for ix = fracElements
        I(1) = 1;
        J(1) = 1;
        V(1) = 1;
        b(1) = newFrac.temperature(ix);
        ind = 2;
        for iy = 2:ny-1
            aL = k / dy(iy) / (yc(iy)-yc(iy-1));
            aR = k / dy(iy) / (yc(iy+1)-yc(iy));
            I(ind) = iy;
            J(ind) = iy-1;
            V(ind) = -aL;
            ind = ind + 1;

            I(ind) = iy;
            J(ind) = iy;
            V(ind) = rho*C/dt + aL + aR;
            ind = ind + 1;

            I(ind) = iy;
            J(ind) = iy+1;
            V(ind) = -aR;
            ind = ind + 1;

            b(iy) = rho*C/dt*oldFrac.rockTemperature.T(ix, iy);
        end

        I(ind) = ny;
        J(ind) = ny;
        V(ind) = 1;
        b(ny) = TInf(ix);

        A = sparse(I, J, V, ny, ny);
        x = A\b';
        
        if covering(ix) == 1
            newFrac.rockTemperature.T(ix,:) = x';
        end
    end
end

