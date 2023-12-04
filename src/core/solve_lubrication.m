function isConvergence = solve_lubrication(startTime, endTime, mesh, reservoir, flowrate, settings, newFrac, oldFrac)
    nIter = 0;
    solError = 1.0;    
    dt = endTime - startTime;
    isConvergence = true;

    
    %% Define channel- and tip-elements
    channel = newFrac.channel;
    tip = newFrac.tip;
    wOldChannel = oldFrac.width(channel);
    wOldTip = oldFrac.width(tip);
    wNewChannel = newFrac.width(channel);
    wNewTip = newFrac.width(tip);
    globalWidth = newFrac.width;

    
    %%
    Q = zeros(1, mesh.n);
    Q(newFrac.perforation) = dt * flowrate / length(newFrac.perforation);
    Qc = Q(channel) ./ mesh.dx(channel);
    Ccc = reservoir.elasticityMatrix(channel, channel);
    Cct = reservoir.elasticityMatrix(channel, tip);
    Ctc = reservoir.elasticityMatrix(tip, channel);
    Ctt = reservoir.elasticityMatrix(tip, tip);
    
    sigmaHc = reservoir.sigmaH(channel);
    sigmaHt = reservoir.sigmaH(tip);
    
    %% Loop by channel width
    sol = zeros(length(channel)+length(tip), 1);
    wIterChannel = wNewChannel;
    while solError > settings.LUBRICATION_TOLERANCE &&...
            nIter < settings.MAX_LUBRICATION_ITERATION
        [Acc, Act, Atc, Att] = generate_mobility_matrix(mesh, globalWidth, channel, tip, reservoir.mup);
        M = assemble_full_matrix(Ccc, Acc, Act, Atc, Att, dt, wIterChannel, settings.MIN_WIDTH);
        rhs = assemble_rhs(Cct, Acc, Act, Atc, Att, dt, wOldChannel, wOldTip, wNewTip, Qc, sigmaHc, sigmaHt, wIterChannel, settings.MIN_WIDTH);
        sol = M\rhs;
        wNewChannel = sol(1:length(channel));
        solError = norm(wNewChannel' - wIterChannel, Inf) / norm(wIterChannel, Inf);
        wIterChannel = wNewChannel';
        nIter = nIter + 1;
        globalWidth(channel) = wIterChannel;
    end
    
    %% Lubrication loop didn't converge
    if nIter >= settings.MAX_LUBRICATION_ITERATION
        fprintf('Lubrication loop didn''t converge');
        isConvergence = false;
        return;
%         error('Lubrication loop didn''t converge');
    end
    
    %% Post-process
    pressureChannel = Ccc * wNewChannel + Cct * wNewTip' + sigmaHc';
    pressureTip = Ctc * wNewChannel + Ctt * wNewTip' + sigmaHt';
    newFrac.set_channel(wNewChannel', pressureChannel');
    newFrac.set_tip_pressure(pressureTip');  
end %solve_lubrication



%% Assemble matrix for solving lubrication equation
function M = assemble_full_matrix(Ccc, Acc, Act, Atc, Att, dt, wIterChannel, wEps)
    penalty = 1e20;
    I = eye(length(Acc));
    wIndicator = penalty * (wIterChannel < wEps);
    Mcc = I - dt * Acc * (Ccc + diag(wIndicator));
    Mct = -dt * Act;
    Mtc = -dt * Atc * Ccc;
    Mtt = -dt * Att;
    M = [Mcc Mct; Mtc Mtt];
end %assemble_full_matrix


%% Assemble RHS in lubrication systems
function rhs = assemble_rhs(Cct, Acc, Act, Atc, Att, dt,...
        wOldChannel, wOldTip, wNewTip, Qc, sigmaHc, sigmaHt, wIterChannel, wEps)
    penalty = 1e20;
    wIndicator = wEps * penalty * (wIterChannel < wEps);
    rhsC = wOldChannel' + dt * Acc * (Cct * wNewTip' + sigmaHc' - wIndicator') + Act * sigmaHt' + Qc';
    rhsT = wOldTip' - wNewTip' + dt * Atc * (Cct * wNewTip' + sigmaHc') + Att * sigmaHt';
    rhs = [rhsC; rhsT];
end %assemble_rhs
