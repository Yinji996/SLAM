function [position_DR,velocity_DR] = Dead_Reckoning(position,gyro_heading_deg)
% parameters
deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor

Dead_reckoning = csvread('Dead_reckoning.csv');
speed_average = (Dead_reckoning(:,2)+Dead_reckoning(:,3)+Dead_reckoning(:,4)+Dead_reckoning(:,5))/4;
heading_deg = gyro_heading_deg;
heading_rad = heading_deg * deg_to_rad;

times = size(gyro_heading_deg,1); %calculate the number of terms

%calculate average velocity
v_N_average = zeros(times,1);
v_E_average = zeros(times,1);
for i = 2:times
    v_N_average(i) = (cos(heading_rad(i)) + cos(heading_rad(i-1)))*speed_average(i)/2;
    v_E_average(i) = (sin(heading_rad(i)) + sin(heading_rad(i-1)))*speed_average(i)/2;
end

%calculate latitude and longitude
latitude_rad = zeros(times,1);
longitude_rad = zeros(times,1);
latitude_rad(1) = position(1,1)*deg_to_rad;
longitude_rad(1) = position(1,2)*deg_to_rad;

for i = 2:times
    [R_N,R_E]= Radii_of_curvature(latitude_rad(i-1));
    latitude_rad(i) = latitude_rad(i-1) + v_N_average(i)*0.5/(R_N+position(i-1,3));
    longitude_rad(i) = longitude_rad(i-1) + v_E_average(i)*0.5/((R_E+position(i-1,3))*cos(latitude_rad(i)));
end

latitude_deg = latitude_rad * rad_to_deg;
longitude_deg = longitude_rad * rad_to_deg;

% calculate the instantaneous DR velocity
v_N_instant = zeros(times,1);
v_E_instant = zeros(times,1);
v_N_instant(1) = speed_average(1)*cos(heading_rad(1));
v_E_instant(1) = speed_average(1)*sin(heading_rad(1));

d = 0.6;
for i = 2:times
    v_N_instant(i) = (2-d)*v_N_average(i) - (1-d)*v_N_instant(i-1);
    v_E_instant(i) = (2-d)*v_E_average(i) - (1-d)*v_E_instant(i-1);
end

position_DR = [latitude_deg,longitude_deg];
velocity_DR = [v_N_instant,v_E_instant];

