classdef OdometryOnlyLocalizationSystem < minislam.localization.ISAMLocalizationSystem
   
    % This class extends the localization system to properly handle vehicle
    % odometry
    
    methods(Access = public)
        
        % Call the base class constructor
        function this = OdometryOnlyLocalizationSystem()
            this = this@minislam.localization.ISAMLocalizationSystem();
        end
        
    end
    
    methods(Access = protected)
        
        % This method should return the relative transformation and the
        % covariance in this relative transformation given a time step length
        % dT. The solution below simply predicts forwards at the speed and
        % ignores steer information. The covariance is also arbitrary.
        function [relativePose, relativePoseCovariance] = computeRelativeVehicleTransformAndCovariance(this, dT)
            
            % The vector this.u contains the most recent control inputs in
            % the form [speed, turn angle rate].
            
            %vDT = dT * this.u(1);
            %relativePose = gtsam.Pose2([dT * this.u(1);dT * this.u(2); 0]);           
            %relativePoseCovariance = diag([0.2 0.1 0.1])*dT^2;
            
            % Compute the relative translation and rotation
            relativeTranslation = gtsam.Point2([dT *this.u(1) * cos(0.5 * this.u(2));dT * this.u(1) * sin(0.5 * this.u(2))]);   
            relativeRotation = dT * this.u(2);
    
            % Assemble the relative transformation.
            relativePose = gtsam.Pose2(relativeRotation, relativeTranslation);
            
            %Calculate system noise covariance matrix
            Qd = zeros(3, 2);
            Qd(1, 1) = dT * cos(0.5 * dT * this.u(2));
            Qd(1, 2) = -0.5 * dT * this.u(1) * sin(0.5 * dT * this.u(2));
            Qd(2, 1) = dT * sin(0.5 * dT * this.u(2));
            Qd(2, 2) = 0.5 * dT * this.u(1) * cos(0.5 * dT * this.u(2));
            Qd(3, 2) = dT;    

            RU = this.uCov;
            relativePoseCovariance = Qd * RU * Qd';
        end
    end
end