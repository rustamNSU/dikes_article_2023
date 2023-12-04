function [Vg, time, dP, rhog, alpha, kg] = gas_filtration_survey(simDir, timesteps)
    %% Load files from folder
    sim_filename = @(i) sprintf('timestep_%05d.mat', i);
    Vg   = zeros(1, length(timesteps));
    time = zeros(1, length(timesteps));
    dP   = zeros(1, length(timesteps));
    rhog = zeros(1, length(timesteps));
    alpha = zeros(1, length(timesteps));
    kg   = zeros(1, length(timesteps));
    time = zeros(1, length(timesteps));
    iter = 1;
    g = 9.81;
    k0 = 5.0e-12;
    j = 3.5;
    mu = 1.5e-5;

    for i = timesteps
        filepath = fullfile(simDir, sim_filename(i));
        load(filepath, 'FractureData');
        if i==100
            a=1;
        end
        
        surveyElement = FractureData.survey(end);
        p1 = FractureData.pressure(surveyElement-1);
        p2 = FractureData.pressure(surveyElement);
        x1 = FractureData.mesh.xc(surveyElement-1);
        x2 = FractureData.mesh.xc(surveyElement);
        dP(iter) = (p2-p1) / (x2-x1);
        time(iter) = FractureData.time;
        rhog(iter) = FractureData.rhogFrac(end);
        alpha(iter) = FractureData.alphaFrac(end);
        kg(iter) = k0 * alpha(iter)^j;
        Vg(iter) = -kg(iter) / mu * (dP(iter) - rhog(iter) * g);
        
        iter = iter+1;
    end
end

