% This script can be used to compare your SLAM system

% Script to set up search paths. 
addpath(genpath('/opt/gtsam/compgx04'))
addpath(genpath(pwd))

% Configure to disable other sensor types
parameters = minislam.event_generators.simulation.Parameters();

% Magic tuning - do not change :)
parameters.laserDetectionRange = 20;

% By setting true / false you can enable different combinations of sensors
parameters.enableOdometry = false;
parameters.enableGPS = true;
parameters.enableLaser = true;

% This changes the period (time between each measurement) of the GPS
parameters.gpsMeasurementPeriod = 500;

% Set up the simulator and the output
eventGenerator = minislam.event_generators.simulation.Simulator(parameters, 'task3');

% Run GTSAM
gtsamSLAMSystem = answers.gtsam.LaserSensor2DSLAMSystem();
gtsamResults = minislam.mainLoop(eventGenerator, gtsamSLAMSystem);

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(gtsamResults.optimizationTimes)
