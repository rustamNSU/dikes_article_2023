clear all
close all
clc

%% ILSA parameters
settingsInput.MIN_WIDTH = 1e-10;
settingsInput.IS_CONSTANT_DENSITY = true;
settingsInput.IS_MAGMA_WEIGHT_ON = true;
settingsInput.IS_CONSTANT_TEMPERATURE = true;
settingsInput.IS_CRYSTALLIZATION_ON = false;
settingsInput.IS_VISCOUS_DISSIPATION_ON = false;
settingsInput.IS_PRESSURE_FORCES_WORK_ON = false;
settingsInput.IS_HOST_TEMPERATURE_ON = false;
settingsInput.FIXED_BOTTOM_TIP = true;
settingsInput.IS_VARIABLE_HOST_CONDUCTIVITY = false;
settingsInput.FRACTURE_TYPE = "KGD";


% "equilibrium", "linear", "distribution", "none"
settingsInput.CRYSTALLIZATION_TYPE = "none";

%% Reservoir parameters
reservoirInput.Young = 20e9;       % [GPa], Young's modulus
reservoirInput.KIc   = 2e6;  % [Pa*m^(1/2)], fracture toughness
reservoirInput.nu    = 0.25;       % [], Poisson's ratio
reservoirInput.rho   = 2000.0;       % [kg/m^3]
reservoirInput.gravity = 10.0;     % [m/s^2]
reservoirInput.h = 35; % [m]


%% Mesh parameters
xmin = -400; % left boundary of domain
xmax = 0;  % right boundary of domain
N = 400;   % number of mesh elements

xi0 = linspace(xmin, xmax, N+1);      % [m], boundary points of elements
xc = (xi0(1:end-1) + xi0(2:end)) / 2; % [m],  middle points
dx = ones(1, N) * (xc(2)-xc(1));

timeList = [0, 1000, 5000, 20000, 100000, 300000];
dtList   = [5, 10, 20, 50, 100]; % [s], time step

reservoirInput.sigmaH = reservoirInput.rho * reservoirInput.gravity * (0 - xc);
xPerforation = -390; % [m], source point
pumpingSchedule = [0 39; 0.05 / reservoirInput.h, 0.0]; % first row: start time, second row: pumping rate
rhoPerforation = 1000; % magma density in chamber
mu = 0.05; % [Pa*s]

filename = 'input/input_simID_5.mat';
save(filename, 'reservoirInput','settingsInput', 'xc', 'dx', 'timeList', 'dtList', 'xPerforation', 'pumpingSchedule', 'rhoPerforation', 'mu');

