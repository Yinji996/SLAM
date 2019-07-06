function [x_init] = Clock_offset_drift()

% define parameters
deg_to_rad = 0.01745329252; 
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

% convert the position to Cartesian ECEF position
latitude = 0;
longitude = 0;
height = 0;
[r_eb_e,v_eb_e] = pv_NED_to_ECEF(latitude*deg_to_rad,longitude*deg_to_rad,height,[0; 0; 0]);

% initialize values
Sat_r_es_e = zeros(8,3);
Sat_v_es_e = zeros(8,3);
position = zeros(1,3);
velocity = zeros(1,3);
clock_offset = 100000;
der_clock_offset = 200;
raj = zeros(8,1);
der_raj = zeros(8,1);

csv_Pseudo_ranges = csvread('Pseudo_ranges.csv');
csv_Pseudo_range_rates = csvread('Pseudo_range_rates.csv');
for epoch = 1:1
    
    % compute the Cartesian ECEF positions of the satellites
    j = csv_Pseudo_ranges(1,2:size(csv_Pseudo_ranges,2));% j:satellite numbers
    for i = 1:size(j,2)
        % sat_r_es_e:ECEF satellite position, sat_v_es_e:ECEF satellite velocity
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(0.5*(epoch-1),j(i));
        Sat_r_es_e(i,:) = sat_r_es_e;%%r_ej
        Sat_v_es_e(i,:) = sat_v_es_e;%%v_ej
    end

    %predict the ranges
    for int = 1:20
        for i = 1:size(j,2)
            raj_tem = 0;
            for iteration = 1:6
                C = [1,omega_ie*raj_tem/c,0;-omega_ie*raj_tem/c,1,0;0,0,1];
                raj_tem = sqrt((C*Sat_r_es_e(i,:)' - r_eb_e)'*(C*Sat_r_es_e(i,:)' - r_eb_e));
                
            end
            raj(i) = raj_tem;
        end

        % compute the line-of-sight unit vector
        u = zeros(8,3);
        for i = 1:size(j,2)
            C = [1,omega_ie*raj(i)/c,0;-omega_ie*raj(i)/c,1,0;0,0,1];
            u(i,:) = (C*Sat_r_es_e(i,:)' - r_eb_e)/raj(i);
        end
        
        % predict the range rates
        omega = [0 -omega_ie 0;omega_ie 0 0;0 0 0];
        for i = 1:size(j,2)
            der_raj_tem = 0;
            C = [1,omega_ie*raj(i)/c,0;-omega_ie*raj(i)/c,1,0;0,0,1];
            for iteration = 1:6
                der_raj_tem = u(i,:)*(C*(Sat_v_es_e(i,:)' + omega*Sat_r_es_e(i,:)')-(v_eb_e+omega*r_eb_e));
            end
            der_raj(i) = der_raj_tem;
        end
        
        %formulate the predicted state vector
        last_col = [1;1;1;1;1;1;1;1];
        x_old = [r_eb_e;clock_offset];
        rho = csv_Pseudo_ranges(epoch+1,2:size(csv_Pseudo_ranges,2))';
        delta_z = rho-raj-clock_offset;
        H_G = [-u last_col];
        x_new = x_old + inv(H_G'*H_G)*H_G'*delta_z;
        
        %update r_ea+ and clock offset
        r_eb_e = x_new(1:3);
        clock_offset = x_new(4);
        
        der_x_old = [v_eb_e;der_clock_offset];
        der_rho = csv_Pseudo_range_rates(epoch+1,2:size(csv_Pseudo_ranges,2))';
        der_delta_z = der_rho-der_raj-der_clock_offset;
        H_G = [-u last_col];
        der_x_new = der_x_old + inv(H_G'*H_G)*H_G'*der_delta_z;
        %update r_ea+ and clock offset
        v_eb_e = der_x_new(1:3);
        der_clock_offset = der_x_new(4);
    end
        
    % outlier detection
    v = (H_G/(H_G'*H_G)*H_G'-eye(8))*delta_z;
    sigma_rho = 5;
    C_v = (eye(8) - H_G/(H_G'*H_G)*H_G')*sigma_rho^2;
    outlier_judge = zeros(8,1);
    Threshold = 6;
    for number = 1:8
        tem = abs(v(number))- sqrt(C_v(number,number))*Threshold;
        if tem>0
        outlier_judge(number) = tem;
        end
    end
           
    % convert the Cartesian position and velocity solution
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    latitude = L_b * rad_to_deg;
    longitude = lambda_b * rad_to_deg;
    position(epoch,:) = [latitude,longitude,h_b];
    velocity(epoch,:) = v_eb_n;
    x_init = [r_eb_e',v_eb_e',clock_offset,der_clock_offset]';
end