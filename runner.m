%% Run simulation of magma dikes propogation from input data
clear all;
close all;
clc;
addpath(genpath(pwd));
warning off;


%% Load input data (pre-saved data)
plotFigure = false;
inputPath        = 'input/input_simID_143001.mat';      % input data
saveFilename     = 'simID_143001';                      % dir path where to save data

load(inputPath, 'reservoirInput', 'settingsInput', 'xc', 'dx',...
    'timeList', 'dtList', 'xPerforation', 'pumpingSchedule', 'rhoPerforation');

load(inputPath, 'magmaInput', 'alphaTemp', 'T0', 'GasDataPath', 'DissolvedGasPath');

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
reservoir = Reservoir(reservoirInput);
clear xc dx reservoirInput;

%% Generate elasticity matrix
reservoir.elasticityMatrix = generate_elasticity_matrix(reservoir.Ep, mesh);

%% Initialize magma state class
state = StateSolver(GasDataPath, DissolvedGasPath, magmaInput, T0);
reservoir.alpha = alphaTemp;
clear magmaInput alphaTemp

oldFrac = FractureElements(mesh);
oldFrac = oldFrac.initialize_zero_fracture(settings, schedule.xPerforation, schedule.startTime);
fracElements = oldFrac.get_fracture_elements();
oldFrac.set_constant_temperature(T0);
oldFrac.pressure = reservoir.sigmaH;

pInChamber = reservoir.sigmaH(oldFrac.perforation);
state.set_u0(pInChamber, T0);

if settings.IS_HOST_TEMPERATURE_ON
    if settings.IS_VARIABLE_HOST_CONDUCTIVITY
        oldFrac.rockTemperature = RockTemperatureVariable(reservoir);
    else
        oldFrac.rockTemperature = RockTemperature(reservoir);
    end
else
    oldFrac.rockTemperature = 0;
end

state.updateEqDensityAndViscosity(oldFrac, settings);
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
save(fullfile(sim_dir, 'state.mat'), 'state');
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
    if plotFigure
        plot_fracture_data(fracture_fig, newFrac, reservoir);
    end
    
    %% Save timestep result
    FractureData = newFrac;
    sim_filename = sprintf('timestep_%05d.mat', timestep);
    save(fullfile(sim_dir, sim_filename), 'FractureData');
end   