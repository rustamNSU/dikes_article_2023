function front_data = parameters_in_front(simDir, timesteps)
    %% Load files from folder
    sim_filename = @(i) sprintf('timestep_%05d.mat', i);
    time = zeros(1, length(timesteps));
    alpha = zeros(1, length(timesteps));
    beta = zeros(1, length(timesteps));
    betaeq = zeros(1, length(timesteps));
    temperature = zeros(1, length(timesteps));
    pressure = zeros(1, length(timesteps));
    rho = zeros(1, length(timesteps));
    mu = zeros(1, length(timesteps));
    ch2od = zeros(1, length(timesteps));
    cco2d = zeros(1, length(timesteps));
    xh2od = zeros(1, length(timesteps));
    xh2o = zeros(1, length(timesteps));
    Mc = zeros(1, length(timesteps));
    Mm = zeros(1, length(timesteps));
    Mg = zeros(1, length(timesteps));
    front = zeros(1, length(timesteps));

    for iter = timesteps
        i = iter + 1;
        filepath = fullfile(simDir, sim_filename(iter));
        load(filepath, 'FractureData');        
        ind = FractureData.survey(end);
        time(i) = FractureData.time;
        alpha(i) = FractureData.alpha(ind);
        beta(i) = FractureData.beta(ind);
        betaeq(i) = FractureData.betaeq(ind);
        temperature(i) = FractureData.temperature(ind);
        pressure(i) = FractureData.pressure(ind);
        rho(i) = FractureData.rho(ind);
        mu(i) = FractureData.mu(ind);
        ch2od(i) = FractureData.ch2o(ind);
        cco2d(i) = FractureData.cco2(ind);
        xh2od(i) = FractureData.xh2od(ind);
        xh2o(i) = FractureData.xh2o(ind);
        Mc(i) = FractureData.Mc(ind);
        Mm(i) = FractureData.Mm(ind);
        Mg(i) = FractureData.Mg(ind);
        front(i) = FractureData.front(end);
    end
    front_data.time = time;
    front_data.alpha = alpha;
    front_data.beta = beta;
    front_data.betaeq = betaeq;
    front_data.temperature = temperature;
    front_data.pressure = pressure;
    front_data.rho = rho;
    front_data.mu = mu;
    front_data.ch2od = ch2od;
    front_data.cco2d = cco2d;
    front_data.xh2od = xh2od;
    front_data.xh2o = xh2o;
    front_data.Mc = Mc;
    front_data.Mm = Mm;
    front_data.Mg = Mg;
    front_data.front = front;
end

