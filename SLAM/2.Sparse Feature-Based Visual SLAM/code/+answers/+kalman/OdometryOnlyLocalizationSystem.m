classdef OdometryOnlyLocalizationSystem < minislam.localization.KalmanFilterLocalizationSystem
    
    % This class extends the localization system to properly handle vehicle
    % odometry
    
    methods(Access = public)
        
        % Call the base class constructor
        function this = OdometryOnlyLocalizationSystem()
            this = this@minislam.localization.KalmanFilterLocalizationSystem();
        end
        
    end
       
    methods(Access = protected)
    
        % This method should return the relative transformation and the
        % covariance in this relative transformation given a time step length
        % dT. The solution below simply predicts forwards at the speed and
        % ignores steer information. The covariance is also arbitrary.
        
        function [xPred, Fd, Qd] = predictMeanJacobianNoise(this, dT)
            
            % The control inputs are in this.u
            %vDT = dT * this.u(1);
            
            % Predict the current position and orientation
            xPred = this.xEst;% The previous mean in the Kalmnan filter is this.xEst
            %xPred(1) = xPred(1) + vDT;
            xPred(1) = xPred(1) + dT * this.u(1) * cos(xPred(3) + 0.5 * dT *this.u(2));
            xPred(2) = xPred(2) + dT * this.u(1) * sin(xPred(3) + 0.5 * dT *this.u(2));
            xPred(3) = xPred(3) + dT * this.u(2);

            % The covariance prediction equation is
            % PPred = Fd * PEst * Fd' + Qd;
            Fd = eye(3);
            Fd(1, 3) = - dT * this.u(1) * sin(this.xEst(3) + 0.5 * dT *this.u(2));
            Fd(2, 3) = dT * this.u(1) * cos(this.xEst(3) + 0.5 * dT *this.u(2));
            
            %Calculate system noise covariance matrix
            Qd = zeros(3, 2);
            Qd(1, 1) = dT * cos(this.xEst(3) + 0.5 * dT *this.u(2));
            Qd(1, 2) = -0.5 * dT * this.u(1) * sin(this.xEst(3) + 0.5 * dT *this.u(2));
            Qd(2, 1) = dT * sin(this.xEst(3) + 0.5 * dT *this.u(2));
            Qd(2, 2) = 0.5 * dT * this.u(1) * cos(this.xEst(3) + 0.5 * dT *this.u(2));
            Qd(3, 2) = dT;
            
            RU = this.uCov;
            Qd = Qd * RU * Qd';
        end
    end
end