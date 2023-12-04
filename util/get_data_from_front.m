function [w, p, rho, mu, T, alpha, beta, cco2, ch2o, xh2o, xh2od, e, xFront, time, mcGas, mcCrystal] = get_data_from_front(simDir, timesteps)
    %% Load files from folder
    sim_filename = @(i) sprintf('timestep_%05d.mat', i);
    w     = zeros(1, length(timesteps));
    p     = zeros(1, length(timesteps));
    rho   = zeros(1, length(timesteps));
    mu    = zeros(1, length(timesteps));
    T     = zeros(1, length(timesteps));
    alpha = zeros(1, length(timesteps));
    beta  = zeros(1, length(timesteps));
    cco2  = zeros(1, length(timesteps));
    ch2o  = zeros(1, length(timesteps));
    xh2o  = zeros(1, length(timesteps));
    xh2od = zeros(1, length(timesteps));
    e     = zeros(1, length(timesteps));
    xFront= zeros(1, length(timesteps));
    time  = zeros(1, length(timesteps));
    mcGas = zeros(1, length(timesteps));
    mcCrystal = zeros(1, length(timesteps));
    
    
    iter = 1;
    for i = timesteps
        filepath = fullfile(simDir, sim_filename(i));
        load(filepath, 'FractureData');        
        frontElement = FractureData.tip(end);
        
        w(iter)     = FractureData.width(frontElement);
        p(iter)     = FractureData.pressure(frontElement);
        rho(iter)   = FractureData.rho(frontElement);
        mu(iter)    = FractureData.mu(frontElement);
        T(iter)     = FractureData.temperature(frontElement);
        alpha(iter) = FractureData.alphaFrac(end);
        beta(iter)  = FractureData.betaFrac(end);
        cco2(iter)  = FractureData.cco2Frac(end);
        ch2o(iter)  = FractureData.ch2oFrac(end);
        xh2o(iter)  = FractureData.xh2oFrac(end);
        xh2od(iter) = FractureData.xh2odFrac(end);
        e(iter)     = FractureData.u(frontElement);
        xFront(iter)= FractureData.front(end);
        time(iter)  = FractureData.time;
        mcGas(iter) = FractureData.MgFrac(end);
        mcCrystal(iter) = FractureData.McFrac(end);
        iter = iter + 1;
    end
end


