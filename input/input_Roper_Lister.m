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


% "equilibrium", "linear", "distribution", "none"
settingsInput.CRYSTALLIZATION_TYPE = "none";

%% Reservoir parameters
reservoirInput.Young = 18.75e9;       % [GPa], Young's modulus
reservoirInput.KIc   = 2.34e8;  % [Pa*m^(1/2)], fracture toughness
reservoirInput.nu    = 0.25;       % [], Poisson's ratio
reservoirInput.rho   = 2700.0;       % [kg/m^3]
reservoirInput.gravity = 10.0;     % [m/s^2]


%% Mesh parameters
xmin = -30000; % left boundary of domain
xmax = 0;  % right boundary of domain
N = 1000;   % number of mesh elements

xi0 = linspace(xmin, xmax, N+1);      % [m], boundary points of elements
xc = (xi0(1:end-1) + xi0(2:end)) / 2; % [m],  middle points
dx = ones(1, N) * (xc(2)-xc(1));

timeList = [0, 3000];
dtList   = [1, 1]; % [s], time step

reservoirInput.sigmaH = reservoirInput.rho * reservoirInput.gravity * (0 - xc);
xPerforation = -29000; % [m], source point
pumpingSchedule = [0 3000; 20.0, 0.0]; % first row: start time, second row: pumping rate
rhoPerforation = 2400; % magma density in chamber
mu = 100; % [Pa*s]

filename = 'input/input_simID_1.mat';
save(filename, 'reservoirInput','settingsInput', 'xc', 'dx', 'timeList', 'dtList', 'xPerforation', 'pumpingSchedule', 'rhoPerforation', 'mu');

