classdef GPSLocalizationSystem < answers.kalman.OdometryOnlyLocalizationSystem
   
    % This class extends the odometry system to process GPS measurements
    % when they are available.
    
    methods(Access = public)
        
        function this = GPSLocalizationSystem()
            this = this@answers.kalman.OdometryOnlyLocalizationSystem();
        end
        
    end
    
    methods(Access = protected)        
        
        % Implement the GPS measurement.
        function handleGPSEvent(this, event)
            
            % You will need to write a Kalman filter update
            
            % The variables this.xPred, this.PPred contain the predicted
            % mean and covariance.
            
            P = minislam.event_generators.simulation.Parameters();
            % store the measurement noise covariance matrix
            RZ = zeros(3:3);
            RZ(1:2,1:2) = P.RGPS;
            RZ(3,3) = 10;
            % Implement a Kalman filter update rule here
            H = eye(3);
            zHat = H * this.xPred;
            data = zeros(3,1);
            data(1:2) = event.data;
            data(3) = this.xPred(3);
            % W:kalman gain(previous K)
            C = this.PPred * H';
            %RZ = diag([0.1; 0.1]);% Covariance of measuring position and orientation
            S  = H * C + RZ;
            W = C / S;
            % Your update should produce a new estimate this.xEst and this.PEst
            this.xEst = this.xPred + W * (data - zHat);
            this.PEst = this.PPred - W * S * W';
            
        end
    end    
end