function statusLubrication = solve_compressibility_lubrication(startTime, endTime, mesh, reservoir, schedule, settings, state, newFrac, oldFrac)
    iter = 0;
    solError = 1.0; 
    rhoError = 1.0;
    dt = endTime - startTime;
    
    if settings.IS_MAGMA_WEIGHT_ON
        isBuoyancy = 1;
    else
        isBuoyancy = 0;
    end

    
    %% Define channel- and tip-elements
    channel = newFrac.channel;
    tip = newFrac.tip;

    %% Set density and viscosity in tip elements from survey
    correctionInTip = tip;
    correctionInTip(1:length(newFrac.leftTip)) = newFrac.survey(1);
    correctionInTip(length(newFrac.leftTip)+1:end) = newFrac.survey(2);

    wOldChannel = oldFrac.width(channel);
    wOldTip = oldFrac.width(tip);
    wNewChannel = newFrac.width(channel);
    wNewTip = newFrac.width(tip);
    
    mOldChannel = oldFrac.rho(channel) .* wOldChannel;
    mOldTip = oldFrac.rho(tip) .* wOldTip;
    rhoNewChannel = newFrac.rho(channel);
    rhoNewTip = newFrac.rho(correctionInTip);
    
    %%
    dMInjected = schedule.injected_mass(startTime, endTime);
    perforation = newFrac.perforation;
    mQ = zeros(1, mesh.n);
    mQCrystal = zeros(1, mesh.n);
    mQ(perforation) = dMInjected / length(perforation) / mesh.dx(perforation);
    Ccc = reservoir.elasticityMatrix(channel, channel);
    Cct = reservoir.elasticityMatrix(channel, tip);
    Ctc = reservoir.elasticityMatrix(tip, channel);
    Ctt = reservoir.elasticityMatrix(tip, tip);
    
    sigmaHc = reservoir.sigmaH(channel);
    sigmaHt = reservoir.sigmaH(tip);

    %% Loop by channel width
    wIterChannel = wNewChannel;
    rhoIterChannel = rhoNewChannel;
    statusLubrication.hasNegativeWidth = false;
    
    while (solError > settings.LUBRICATION_TOLERANCE &&...
           iter < settings.MAX_LUBRICATION_ITERATION)
        mQCrystal(perforation) = newFrac.rhoc(perforation) * mQ(perforation) / newFrac.rho(perforation);
        magmaWeightGlobal = isBuoyancy * reservoir.gravity * newFrac.rho .* mesh.xc;
        [Acc, Act, Atc, Att] = generate_compressibility_mobility_matrix(mesh, newFrac);
        
        %% Solving mass law equation
        M = assemble_full_matrix(Ccc, Acc, Act, Atc, Att, dt, wIterChannel, rhoIterChannel, settings.MIN_WIDTH);
        rhs = assemble_rhs(Cct, Acc, Act, Atc, Att, dt,...
            mOldChannel, mOldTip, wNewTip, rhoNewTip, mQ(channel), ...
            sigmaHc, sigmaHt, wIterChannel, settings.MIN_WIDTH,...
            magmaWeightGlobal(channel), magmaWeightGlobal(tip));
        sol = M \ rhs;
        
        widthRelaxation = settings.get_lubrication_relaxation_coef(iter);
        wNewChannel = sol(1:length(channel))';
        if any(wNewChannel < 0) || any(isnan(wNewChannel))
            statusLubrication.hasNegativeWidth = true;
            fprintf(2, 'Negative width in solve_compressibility_lubrication\n');
            break;
        end
        wNewChannel = (1.0 - widthRelaxation) * wIterChannel + widthRelaxation * wNewChannel;
        widthError = norm(wNewChannel - wIterChannel, Inf) / norm(wIterChannel, Inf);
        wIterChannel = wNewChannel;
        iter = iter + 1;
        newFrac.width(channel) = wNewChannel;
        
        %% Find mass rate and flux velocity
        find_mass_rate(newFrac, oldFrac, mesh, mQ, dt);
        globalTIter = newFrac.temperature;
        if not(settings.IS_CONSTANT_TEMPERATURE) 
            solve_energy_conservation(mesh, newFrac, oldFrac,...
            dt, mQ, mQCrystal, reservoir, settings, state);
            newFrac.temperature = widthRelaxation * newFrac.temperature + (1.0 - widthRelaxation) * globalTIter;
        end

        %% Update pressure
        pressureChannel = (Ccc * wNewChannel' + Cct * wNewTip')' + sigmaHc;
        pressureChannel(pressureChannel <= 0) = 0;
        % if any(pressureChannel < 0) || any(pressureChannel > 1.5e9)
        %     statusLubrication.hasNegativeWidth = true;
        %     fprintf(2, 'Pressure in channel not in [0.0, 1.5e9] Pa\n');
        %     break;
        % end
        newFrac.pressure(channel) = pressureChannel;
        newFrac.pressure(tip) = newFrac.pressure(correctionInTip);
        
        %% Update magma state in fracture (set neighbor survey element data in tip)
        rhoError = 0;
        TError = 0;
        if not(settings.IS_CONSTANT_DENSITY) 
            state.updateDensityAndViscosity(newFrac, oldFrac, startTime, endTime, mQ, settings);
            rhoNewChannel  = newFrac.rho(channel);
            rhoNewTip      = newFrac.rho(tip);
            rhoError       = norm(rhoNewChannel - rhoIterChannel, Inf) / norm(rhoIterChannel, Inf);
            rhoIterChannel = rhoNewChannel;
            TError = norm(newFrac.temperature(channel) - globalTIter(channel), Inf) / norm(newFrac.temperature(channel), Inf);
        end
        solError = max([widthError, rhoError, TError]);
    end
    
    statusLubrication.iterations = iter;
    statusLubrication.isNanSolution = isnan(solError);
    
    %% Lubrication loop didn't converge
    if iter >= settings.MAX_LUBRICATION_ITERATION || statusLubrication.isNanSolution || statusLubrication.hasNegativeWidth
        statusLubrication.isConvergence = false;
        return;
    end
     
    statusLubrication.isConvergence = true;
    statusLubrication.solError = solError;
    statusLubrication.rhoError = rhoError;
    statusLubrication.widthError = widthError;
    statusLubrication.tempError = TError;
end %solve_lubrication



%% Assemble matrix for solving lubrication equation
function M = assemble_full_matrix(Ccc, Acc, Act, Atc, Att, dt, wIterChannel, rhoIterChannel, wEps)
    penalty = 1e20;
    rhoI = diag(rhoIterChannel);
    wIndicator = penalty * (wIterChannel < wEps);
    Mcc = rhoI - dt * Acc * (Ccc + diag(wIndicator));
    Mct = -dt * Act;
    Mtc = -dt * Atc * Ccc;
    Mtt = -dt * Att;
    M = [Mcc Mct; Mtc Mtt];
end %assemble_full_matrix


%% Assemble RHS in lubrication systems
function rhs = assemble_rhs(Cct, Acc, Act, Atc, Att, dt,...
        mOldChannel, mOldTip, wNewTip, rhoNewTip,...
        mQc, sigmaHc, sigmaHt, wIterChannel, wEps, magmaWeightChannel, magmaWeightTip)
    penalty = 1e20;
    wIndicator = wEps * penalty * (wIterChannel < wEps);
    rhsC = mOldChannel' + dt * Acc * (Cct * wNewTip' + magmaWeightChannel' + sigmaHc' - wIndicator') +...
        Act * (sigmaHt' + magmaWeightTip') + mQc';
    rhsT = mOldTip' - rhoNewTip'.*wNewTip' + dt * Atc * (Cct * wNewTip' + sigmaHc' + magmaWeightChannel') +...
        Att * (sigmaHt' + magmaWeightTip');
    rhs = [rhsC; rhsT];
end %assemble_rhs


%% State equation
function rho = get_density_linear(pressure, Cf, rho0)
    rho = rho0 + Cf*pressure;
    rho(rho<rho0) = rho0;
    rhoMax = 2300;
    rho(rho>rhoMax) = rhoMax;
end
    