%% Run simulation of magma dikes propogation from input data
clear all;
close all;
clc;
addpath(genpath(pwd));
warning off;


%% Load input data (pre-saved data)
plotFigure = false;
inputPath        = 'input/input_simID_1.mat';      % input data
saveFilename     = 'simID_1';                      % dir path where to save data

load(inputPath, 'reservoirInput', 'settingsInput', 'xc', 'dx',...
    'timeList', 'dtList', 'xPerforation', 'pumpingSchedule', 'rhoPerforation');

if settingsInput.IS_CONSTANT_DENSITY
    load(inputPath, 'mu');
else
    load(inputPath, 'magmaInput', 'alphaTemp', 'T0', 'GasDataPath', 'DissolvedGasPath');
end

%% Save simulation data
sim_dir = fullfile(pwd, 'Simulations', saveFilename);
if not(isfolder(sim_dir))
    mkdir(sim_dir);
end
%% Initialize simulation settings
settings = IlsaSettings(settingsInput);
settings.simlog = fopen(fullfile(sim_dir, 'simlog.txt'), 'wt');
fprintf(settings.simlog, 'inputPath = %s\nsaveFilename = %s\n\n', inputPath, saveFilename);

%% Magma injection from chamber
schedule = Schedule(timeList, dtList, pumpingSchedule, xPerforation, rhoPerforation);
time = schedule.startTime;
clear timeList dtList pumpingSchedule xPerforation rhoPerforation;

%% Initialize mesh and reservoir parameters
mesh = Mesh(xc, dx);
reservoir = Reservoir(reservoirInput, settings);
clear xc dx reservoirInput;

%% Generate elasticity matrix
reservoir.elasticityMatrix = generate_elasticity_matrix(reservoir.Ep, mesh);

oldFrac = FractureElements(mesh);
oldFrac = oldFrac.initialize_zero_fracture(settings, schedule.xPerforation, schedule.startTime);
fracElements = oldFrac.get_fracture_elements();
oldFrac.pressure = reservoir.sigmaH;

state = 0; % dummy value
if settingsInput.IS_CONSTANT_DENSITY == false
    %% Initialize magma state class
    state = StateSolver(GasDataPath, DissolvedGasPath, magmaInput, T0);
    reservoir.alpha = alphaTemp;
    clear magmaInput alphaTemp
    oldFrac.set_constant_temperature(T0);

    pInChamber = reservoir.sigmaH(oldFrac.perforation);
    state.set_u0(pInChamber, T0);
end

if settings.IS_HOST_TEMPERATURE_ON
    if settings.IS_VARIABLE_HOST_CONDUCTIVITY
        oldFrac.rockTemperature = RockTemperatureVariable(reservoir);
    else
        oldFrac.rockTemperature = RockTemperature(reservoir);
    end
else
    oldFrac.rockTemperature = 0;
end

if settingsInput.IS_CONSTANT_DENSITY == true
    oldFrac.mu(:) = mu;
    oldFrac.rho(:) = schedule.rhoPerforation;
    settings.lubrication_relaxation_iter = [20 40 70 100 200];
    settings.lubrication_relaxation_coef = [1.0 0.7 0.5 0.3 0.1];
    settings.MAX_LUBRICATION_ITERATION = 40;
else
    state.updateEqDensityAndViscosity(oldFrac, settings);
    save(fullfile(sim_dir, 'state.mat'), 'state');
end
oldFrac.crystallization = 0;
newFrac = copy(oldFrac);

%% Create figure window for plotting fracture after nonlinear loop
if plotFigure
    fracture_fig = figure('Name','Fracture parameters', 'Color','white',...
        'Units','normalized', 'Position',[0.1, 0.1, 0.6, 0.6]);
end

FractureData = newFrac;
timestep = 0;
sim_filename = sprintf('timestep_%05d.mat', timestep);
save(fullfile(sim_dir, sim_filename), 'FractureData');
save(fullfile(sim_dir, 'reservoir.mat'), 'reservoir');
save(fullfile(sim_dir, 'settings.mat'), 'settings');
save(fullfile(sim_dir, 'schedule.mat'), 'schedule');

%% Start fracture simulation
while time < schedule.endTime - 0.1
    startTime = time;
    time = startTime + schedule.get_dt(startTime);
    timestep = timestep + 1;

    compute_time_step(startTime, time,...
        mesh, reservoir, schedule, settings, state, newFrac, oldFrac, 0);
    
    %% Plot timestep result
    if plotFigure && not(settingsInput.IS_CONSTANT_DENSITY)
        plot_fracture_data(fracture_fig, newFrac, reservoir);
    end
    
    %% Save timestep result
    FractureData = newFrac;
    sim_filename = sprintf('timestep_%05d.mat', timestep);
    save(fullfile(sim_dir, sim_filename), 'FractureData');
end   