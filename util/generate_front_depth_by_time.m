function [front, time] = generate_front_depth_by_time(simDir, timesteps)
    %% Load files from folder
    sim_filename = @(i) sprintf('timestep_%05d.mat', i);
    front = zeros(1, length(timesteps));
    time = zeros(1, length(timesteps));
    iter = 1;
    for i = timesteps
        filepath = fullfile(simDir, sim_filename(i));
        load(filepath, 'FractureData');

        front(iter) = FractureData.front(2);
        time(iter) = FractureData.time;
        iter = iter+1;
    end
end

