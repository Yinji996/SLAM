classdef LaserSensor2DSLAMSystem < answers.gtsam.GPSLocalizationSystem
    
    methods(Access = public)
        
        function this = LaserSensor2DSLAMSystem()
            this = this@answers.gtsam.GPSLocalizationSystem();
        end
    end
    
    methods(Access = protected)
        
        
        function handleLaserEvent(this, event)
            
            % Iterate through each landmark in turn
            numLandmarks = length(event.landmarkIDs);
            
            for l = 1 : numLandmarks
                this.processLaserObservation(event.landmarkIDs(l), event.data(:, l), event.covariance);
            end
            
        end
        
        function processLaserObservation(this, id, z, R)
            
            % The laser scan data works with 3D landmarks (x, y, z) but,
            % unfortunately, the mapping system only works in 2D.
            % Therefore, we have to convert from 3D to 2D to get the range
            % correct.
            S = [1 0 -z(1)*sin(z(3));0 1 0];
            R = S * R * S';% measurement noise covariance                       
            z = [z(1)*cos(z(3));z(2)];% observation

            % This is the way to find a landmark if we have it registered
            % aready
            if (this.landmarkIDKeyStore.contains(id) == true)
                landmarkKey = this.landmarkIDKeyStore.get(id);
            else
                % Create the new key and insert it into the values map
                landmarkKey = this.getLandmarkKey(id);
                this.landmarkIDKeyStore.insert(id, landmarkKey);
                
                % Work out the intial value. The code will be of the form:
                % initialValue = gtsam.Point2(intialX, initalY);
                initialValue = gtsam.Point2(this.currentVehiclePose.x()+z(1)*cos(z(2)+ this.currentVehiclePose.theta()), this.currentVehiclePose.y()+z(1)*sin(z(2)+ this.currentVehiclePose.theta()));

                this.addNewVariable(landmarkKey,initialValue);
            end

            % Constuct the laser observation model.
            laserNoiseModel = gtsam.noiseModel.Gaussian.Covariance(R);

            % Construct the observation factor
            % Add bearing/range measurement factors
            landmarkObservationFactor = gtsam.BearingRangeFactor2D(this.currentVehiclePoseKey, landmarkKey, gtsam.Rot2(z(2)), z(1), laserNoiseModel);
            % Add the new observation factor to the graph
            this.addNewFactor(landmarkObservationFactor);
        end
    end
end
