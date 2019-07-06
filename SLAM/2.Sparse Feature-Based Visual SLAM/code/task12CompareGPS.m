% This script runs the GPS
import gtsam.*

% Script to set up search paths. 
addpath(genpath('/opt/gtsam/compgx04'))
addpath(genpath(pwd))

% Configure to disable other sensor types
parameters = minislam.event_generators.simulation.Parameters();
parameters.enableGPS = true;
parameters.enableLaser = false;

% Set up the simulator and the output
simulator = minislam.event_generators.simulation.Simulator(parameters);

kalmanFilterGPSSystem = answers.kalman.GPSLocalizationSystem();
kalmanResults = minislam.mainLoop(simulator, kalmanFilterGPSSystem);

% Run the two simulations
gtsamGPSSystem = answers.gtsam.GPSLocalizationSystem();
gtsamResults = minislam.mainLoop(simulator, gtsamGPSSystem);

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(kalmanResults.optimizationTimes)
hold on
plot(gtsamResults.optimizationTimes)
legend('kalmanResults','gtsamResults')
xlabel('Step number');
ylabel('Time required to run the optimiser(s)');

% Plot convariance 
figure;
KalmanCovariance = sum(abs(kalmanResults.vehicleCovarianceHistory).^2,1).^(1/2);
plot(KalmanCovariance);
hold on
GTSAMCovariance = sum(abs(gtsamResults.vehicleCovarianceHistory(1:2,:)).^2,1).^(1/2); 
plot(GTSAMCovariance);
legend('kalmanResults','gtsamResults')
xlabel('Step number');
ylabel('Norm of convariance');
