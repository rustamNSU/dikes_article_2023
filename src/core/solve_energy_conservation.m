function solve_energy_conservation(mesh, newFrac, oldFrac, dt, mQ, mQCrystal, reservoir, settings, state)
    alpha = reservoir.alpha;  % J/(m*K)
    Lm    = state.Lc;
    
    fracElements = newFrac.get_fracture_elements();
    
    %% Turn-off heat conduction in tip elements
    covering = zeros(1, mesh.n);
    covering(fracElements) = 1;
    covering(newFrac.tip) = 0;    
    Qinj = mQ .* state.uInChamber;
    
    % delay crystallization while Q > 0
    if newFrac.time > -1
        crystallization_delay = 1;
    else
        crystallization_delay = 0;
    end
    
    S = generate_energy_mobility_matrix(mesh, newFrac, dt);
    Inew = generate_mass_energy_matrix(newFrac, fracElements);
    M = Inew + S;
    
    Iold = generate_mass_energy_matrix(oldFrac, fracElements);
    QCrystallization = crystallization_delay * settings.get_crystallization_coef() * Lm * crystal_mass_transfer(newFrac, oldFrac, fracElements, mesh, dt, mQCrystal);
    W_pressure = crystallization_delay * settings.get_pressure_work_coef() * get_pressure_work(newFrac, oldFrac, fracElements, mesh, dt, reservoir);
    W_vis = crystallization_delay * settings.get_viscous_dissipation_coef() * get_viscousity_dissipation(newFrac, fracElements, mesh, dt);
    rhs_1 = Iold * oldFrac.u(fracElements)' + Qinj(fracElements)' + QCrystallization' + W_pressure' + W_vis';
    
    % error = 1;
    % while (error > 1e-8)

%     T = newFrac.temperature(fracElements);
%     Ql = 2 * alpha * dt * covering(fracElements) .* (reservoir.temperatureH(fracElements) - T);
%     rhs = rhs_1 + Ql';
    if settings.IS_HOST_TEMPERATURE_ON
        newFrac.rockTemperature.Tm = newFrac.temperature;
        newFrac.rockTemperature.channel = newFrac.channel;
%         newFrac.rockTemperature.channel = fracElements;
        newFrac.rockTemperature.solve_heat_equation(oldFrac.rockTemperature, dt);
        ql = newFrac.rockTemperature.get_heat_flow();
        Ql = covering(fracElements)' .* (-2 * dt * ql(fracElements));
%         Ql = -2 * dt * ql(fracElements);
        rhs = rhs_1 + Ql;
    else
        T = newFrac.temperature(fracElements);
        Ql = 2 * alpha * dt * covering(fracElements) .* (reservoir.temperatureH(fracElements) - T);
        Ql(Ql>0) = 0.0;
        rhs = rhs_1 + Ql';
    end
    sol = M \ rhs;
    newFrac.u(fracElements) = sol';
    if any(isnan(sol))
        error('Stop calculation: nan solution in energy equation\n');
    end

    correctionInTip = newFrac.tip;
    correctionInTip(1:length(newFrac.leftTip)) = newFrac.survey(1);
    correctionInTip(length(newFrac.leftTip)+1:end) = newFrac.survey(2);
    newFrac.u(newFrac.tip)         = newFrac.u(correctionInTip);

    T = state.findTemperature(newFrac);
    newFrac.temperature(fracElements) = T;
end


function S = generate_energy_mobility_matrix(mesh, frac, dt)
    %% Fill global mobility matrix
    A = zeros(mesh.n, mesh.n);
    fracElements = frac.get_fracture_elements();
    for i = 1:length(fracElements)
        element = fracElements(i);
        
        mql = frac.mass_rate(element) * dt;
        mqr = frac.mass_rate(element+1) * dt;
        dxc = mesh.dx(element);
        
        %% Define energy flux
        if mql > 0
            iL = element-1;
        else
            iL = element;
        end
        
        if mqr > 0
            iR = element;
        else
            iR = element+1;
        end
        
        %% Left tip element (boundary)
        if i==1 
            A(element, iR) = A(element, iR) + mqr / dxc;
            continue
        end
         
        %% Right tip element (boundary)
        if i==length(fracElements)
            A(element, iL) = A(element, iL) - mql / dxc;
            continue
        end
              
        %% Inner elements
        A(element, iL) = A(element, iL) - mql / dxc;
        A(element, iR) = A(element, iR) + mqr / dxc;
    end %Fill global mobility matrix
    
    S = sparse(A(fracElements, fracElements));
end


function I = generate_mass_energy_matrix(frac, fracElements)
    n = length(fracElements);
    I = 1:n;
    J = 1:n;
    V = zeros(1, n);
    for i = 1:n
        ind = fracElements(i);
        V(i) = frac.rho(ind) * frac.width(ind);
    end
    I = sparse(I, J, V, n, n);
end


function crystal_mass_transfer = crystal_mass_transfer(newFrac, oldFrac, fracElements, mesh, dt, mQCrystal)
    crystal_mass_transfer = zeros(1, length(fracElements));
    for i = 1:length(fracElements)
        ind = fracElements(i);
        dx = mesh.dx(ind);
        mNew = newFrac.width(ind) * newFrac.rhoc(ind) * dx;
        mOld = oldFrac.width(ind) * oldFrac.rhoc(ind) * dx;
        vLeft = newFrac.flux_velocity(ind);
        vRight = newFrac.flux_velocity(ind+1);
        if vLeft > 0
            qLeft = vLeft * newFrac.rhoc(ind-1);
        else
            qLeft = vLeft * newFrac.rhoc(ind);
        end

        if vRight > 0
            qRight = vRight * newFrac.rhoc(ind);
        else
            qRight = vRight * newFrac.rhoc(ind+1);
        end
%         qLeft = vLeft * 0.5 * (newFrac.rhoc(ind-1) + newFrac.rhoc(ind));
%         qRight = vRight * 0.5 * (newFrac.rhoc(ind) + newFrac.rhoc(ind+1));
        
        crystal_mass_transfer(i) = (mNew - mOld - (qLeft - qRight) * dt) / dx - mQCrystal(ind);
    end
end


function W = get_pressure_work(newFrac, oldFrac, fracElements, mesh, dt, reservoir)
    W = zeros(1, length(fracElements));
    for i = 1:length(fracElements)
        ind = fracElements(i);
        dx = mesh.dx(ind);
        qLeft = newFrac.flux_velocity(ind);
        qRight = newFrac.flux_velocity(ind+1);
        wNew = newFrac.width(ind);
        wOld = oldFrac.width(ind);
        W(i) = -newFrac.pressure(ind) * (dt * (qRight - qLeft) / dx + wNew - wOld);
    end
end


function W_vis = get_viscousity_dissipation(newFrac, fracElements, mesh, dt)
    W_vis = zeros(1, length(fracElements));
    g = 9.81;
    pf = newFrac.pressure + g * newFrac.rho .* newFrac.mesh.xc;
    pf = pf(fracElements);
    dpdx = zeros(1, length(pf)-1);
    for i = 1:length(dpdx)
        ind = fracElements(i);
        dpdx(i) = (pf(i+1) - pf(i)) / newFrac.mesh.dx(ind);
    end
    dpdx = [dpdx(1), dpdx, dpdx(end)];
    for i = 1:length(fracElements)
        ind = fracElements(i);
        dx = mesh.dx(ind);
        qLeft = newFrac.flux_velocity(ind);
        qRight = newFrac.flux_velocity(ind+1);
        dpdxLeft = dpdx(i);
        dpdxRight = dpdx(i+1);
%         q2 = (0.5 * (qLeft + qRight))^2;
%         mu = newFrac.mu(ind);
%         w3 = newFrac.width(ind)^3;
%         W_vis(i) = dt * 12 * mu * q2 / w3;
        W_vis(i) = -dt * (qLeft * dpdxLeft + qRight * dpdxRight) / 2;
    end
end