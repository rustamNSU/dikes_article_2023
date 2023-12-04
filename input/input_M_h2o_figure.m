clear all
close all
clc

% simIDs = {143011, 143021, 143031, 143041, 143051};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% }; % first row: start time, second row: pumping rate
% DissolvedGasPath = 'statePinatubo/pdpath0.1.mat';   % state data for dissolved gas
% timeList = [0, 50000, 150000, 350000, 750000, 1550000];
% dtList   = [50, 100, 200, 400, 800]; % [s], time step


% simIDs = {143012, 143022, 143032, 143042, 143052};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.2.mat';
% timeList = [0, 50000, 150000, 350000, 750000];
% dtList   = [50, 100, 200, 400];

% simIDs = {143013, 143023, 143033, 143043, 143053};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.3.mat';
% timeList = [0, 50000, 150000, 350000];
% dtList   = [50, 100, 200];

% simIDs = {143014, 143024, 143034, 143044, 143054};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.4.mat';
% timeList = [0, 50000, 150000];
% dtList   = [50, 100];

% simIDs = {143015, 143025, 143035, 143045, 143055};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.5.mat';
% timeList = [0, 50000, 150000];
% dtList   = [50, 100];

% simIDs = {143016, 143026, 143036, 143046, 143056};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.6.mat';
% timeList = [0, 50000];
% dtList   = [50];

% simIDs = {143017, 143027, 143037, 143047, 143057};
% pumpingSchedules = {
%     [0 10000; 2, 0.0],...
%     [0 10000; 3, 0.0],...
%     [0 10000; 4, 0.0],...
%     [0 10000; 5, 0.0],...
%     [0 10000; 6, 0.0]
% };
% DissolvedGasPath = 'statePinatubo/pdpath0.7.mat';
% timeList = [0, 50000];
% dtList   = [50];

simIDs = {143019, 143029, 143039, 143049, 143059};
pumpingSchedules = {
    [0 10000; 2, 0.0],...
    [0 10000; 3, 0.0],...
    [0 10000; 4, 0.0],...
    [0 10000; 5, 0.0],...
    [0 10000; 6, 0.0]
};
DissolvedGasPath = 'statePinatubo/pdpath0.9.mat';
timeList = [0, 50000];
dtList   = [50];

for i=1:length(simIDs)
    simID = simIDs{i};
    pumpingSchedule = pumpingSchedules{i};
    generate_input(simID, DissolvedGasPath, pumpingSchedule, timeList, dtList);
end

function generate_input(simID, DissolvedGasPath, pumpingSchedule, timeList, dtList)
    %% ILSA parameters
    settingsInput.MIN_WIDTH = 1e-10;
    settingsInput.IS_CONSTANT_DENSITY = false;
    settingsInput.IS_MAGMA_WEIGHT_ON = true;
    settingsInput.IS_CONSTANT_TEMPERATURE = false;
    settingsInput.IS_CRYSTALLIZATION_ON = true;
    settingsInput.IS_VISCOUS_DISSIPATION_ON = true;
    settingsInput.IS_PRESSURE_FORCES_WORK_ON = false;
    settingsInput.IS_HOST_TEMPERATURE_ON = true;
    settingsInput.FIXED_BOTTOM_TIP = true;
    settingsInput.IS_VARIABLE_HOST_CONDUCTIVITY = true;


    % "equilibrium", "linear", "distribution"
    settingsInput.CRYSTALLIZATION_TYPE = "equilibrium";



    %% Reservoir parameters
    reservoirInput.Young = 15e9; % [GPa], Young's modulus
    reservoirInput.KIc   = 1e6;  % [Pa*m^(1/2)], fracture toughness
    reservoirInput.nu    = 0.25; % [], Poisson's ratio
    reservoirInput.rho = 2700.0; % [kg/m^3]
    reservoirInput.gravity = 9.81;        % [m/s^2]
    reservoirInput.cp = 1200; % [J/(kg*K*s)]
    reservoirInput.k  = 1.9; % [J/(K*m*s) = W/(K*m)], thermal conductivity 
    reservoirInput.Ly = 2.0;
    reservoirInput.Ny = 30;


    %% Magma parameters
    magmaInput.Cv_c = 1200;   % [J/(kg*K)], crystal specific heat
    magmaInput.Cv_m = 1200;   % [J/(kg*K)], melt specific heat
    magmaInput.Lc = 350000;   % [J/kg], latent heat due crystallization
    magmaInput.tau = 10000;     % [s], relaxation time in the crystallization law
    alphaTemp = 10;  % [J/(m^2*K)]


    %% Mesh parameters
    xmin = -33000; % left boundary of domain
    xmax = 200;  % right boundary of domain
    N = 400;   % number of mesh elements

    xi0 = linspace(xmin, xmax, N+1);      % [m], boundary points of elements
    xc = (xi0(1:end-1) + xi0(2:end)) / 2; % [m],  middle points
    dx = ones(1, N) * (xc(2)-xc(1));
    % timeList = [0,   20000, 60000, 200000, 400000, 800000, 1600000];
    % dtList   = [50,  100,  200, 400, 800, 1600]; % [s], time step

    reservoirInput.sigmaH = reservoirPressure(xc);
    xPerforation = -32500; % [m], source point
    rhoPerforation = 2500;


    T0 = 850 + 273.15;
    % T0 = 900 + 273.15;
    dTdz = 40e-3;
    omega = 0.1;
    D1 = 10000;
    Tsurf = 0; 
    temperatureH = reservoirTemperature(xc, dTdz, omega, T0-273.15-Tsurf, D1, xPerforation) + 273.15 + Tsurf;
    temperatureH(temperatureH > T0) = T0;
    reservoirInput.temperatureH = temperatureH;


    GasDataPath      = 'statePinatubo/propH2OCO2.mat'; % from VESICAL package

    filename = "input/input_simID_" + int2str(simID) +".mat";
    save(filename, 'reservoirInput', 'xc', 'dx', 'timeList', 'dtList', 'xPerforation', 'pumpingSchedule', 'rhoPerforation');
    save(filename, 'settingsInput', 'magmaInput', 'alphaTemp', 'T0', 'GasDataPath', 'DissolvedGasPath', '-append');
end



function T = calculate_temperature_by_depth(T0, x0, x)
    Ts = 20 + 273.15; % surface temperature (x = 0)
    xs = 0;
    T = T0 + (x-x0) * (Ts - T0) / (xs - x0);    
end

function T=reservoirTemperature(depth, dTdz, omega, Td, D1, chamberDepth)
    z = -depth;
    b=2e-3;
    H=-chamberDepth;  
    t2 = (D1 ^ 2);
    t4 = 1 / D1;
    t6 = exp(-t4 * H);
    t15 = log(1 / (Td * b + 1));
    t19 = exp(-t4 * z);
    t22 = z .^ 2;
    t24 = 2 * t2;
    t39 = H ^ 2;
    t46 = exp(1 ./ (2 * D1 * H + omega * t39 + 2 * t6 * t2 - t24) * ...
        (2 * t6 * dTdz * t2 * b * z - 2 * t19 * t2 * (H * b * dTdz + t15) + ...
        t15 * (-2 * D1 * z - t22 * omega + t24) + (H * omega * z + t24) .* (H - z) * b * dTdz));
    T = 1 / b * (t46 - 1);
end


% get density of reservoir (z is a depth in [m] with minus sign)
function rho = reservoirDensity(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    d = 1724.7;
    alpha = 0.001918;
    rho = rho1 - alpha * z - (rho1 - rho0) * exp(z / d);
end

% get lithostatic pressure of reservoir (z is a depth in [m] with minus sign)
function p = reservoirPressure(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    g = 9.81;
    d = 1724.7;
    alpha = 0.001918;
    p = g * (-rho1*z + 0.5*alpha*z.*z - d*(rho1-rho0)*(1.0-exp(z/d)));
end