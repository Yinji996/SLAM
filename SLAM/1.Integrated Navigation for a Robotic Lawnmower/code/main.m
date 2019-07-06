clear all;
% obtain the initial value of position, velocity, offset and drift
x_init = Clock_offset_drift(); 
% use extended Kalman Filter to calculate GNSS position and velocity
[outlier_judge,position,velocity] = GNSS_extended_Kalman_filter(x_init);
% Gyro-Magnetometer corrected
gyro_heading_deg = Gyro_heading();
% calculate Dead Reckoning position and velocity
[position_DR,velocity_DR] = Dead_Reckoning(position,gyro_heading_deg);
% DR/GNSS integration
results = Itegration(position,velocity,position_DR,velocity_DR,gyro_heading_deg);
% create motion profile file
csvwrite('Output_Profile.csv',results);