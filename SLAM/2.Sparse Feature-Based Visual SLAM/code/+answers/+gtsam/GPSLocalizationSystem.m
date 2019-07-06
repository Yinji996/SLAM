classdef GPSLocalizationSystem < answers.gtsam.OdometryOnlyLocalizationSystem
   
    % This class extends the odometry system to process GPS measurements
    % when they are available.
    
    methods(Access = public)
        
        function this = GPSLocalizationSystem()
            this = this@answers.gtsam.OdometryOnlyLocalizationSystem();
        end
        
    end
    
    methods(Access = protected)        
        
        % Implement the GPS measurement.
        function handleGPSEvent(this, event)
            
            % You will need to create a suitable factor and noise model
            % insert it into the graph. You will need:
            
            % Insert the initial estimate for this new pose key using the relative
            % transformation we computed above.
            
                      
            % Specify the measurement noise model, e.g.,
            P = minislam.event_generators.simulation.Parameters();
            % store the measurement noise covariance matrix
            gpsCovariance = zeros(3);
            gpsCovariance(1:2,1:2) = P.RGPS;
            gpsCovariance(3,3) = 10;
            gpsNoise = gtsam.noiseModel.Gaussian.Covariance(gpsCovariance);
            
            priorMean = gtsam.Pose2(this.currentVehiclePose.theta(), gtsam.Point2(event.data));
            gpsObservationFactor = gtsam.PriorFactorPose2(this.currentVehiclePoseKey, priorMean, gpsNoise);
            
            % Construct the suitable factor
            %gpsObservationFactor = ...
            
            % Insert it into the graph. Use this method:
            this.addNewFactor(gpsObservationFactor);
        end
    end    
end