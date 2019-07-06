function [outlier_judge,position,velocity] = GNSS_extended_Kalman_filter(x_init)
% parameters
deg_to_rad = 0.01745329252; 
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
c = 299792458; % Speed of light in m/s

% initialize x and P
x_est = x_init;
P_matrix = [eye(3)*10^2,zeros(3,5);
            zeros(3,3),eye(3)*0.1^2,zeros(3,2);
            zeros(1,6),10^2,0;
            zeros(1,7),0.1^2];

% compute the transition matrix
tau_s = 0.5;
I3 = eye(3);
zero3_1 = zeros(3,1);
zero1_3 = zeros(1,3);
phi = [I3,tau_s*I3,zero3_1,zero3_1;
       zeros(3,3),I3,zero3_1,zero3_1;
       zero1_3,zero1_3,1,tau_s;
       zero1_3,zero1_3,0,1];

% compute the system noise convariance matrix
S_a = 5;
S_c_phi = 0.01;
S_c_f = 0.04;
Q = [S_a*tau_s^3*I3/3,S_a*tau_s^2*I3/2,zero3_1,zero3_1;
     S_a*tau_s^2*I3/2,S_a*tau_s*I3,zero3_1,zero3_1;
     zero1_3,zero1_3,S_c_phi*tau_s+S_c_f*tau_s^3/3,S_c_f*tau_s^2/2;
     zero1_3,zero1_3,S_c_f*tau_s^2/2,S_c_f*tau_s];

position = zeros(851,3);
velocity = zeros(851,3);

outlier_judge = zeros(851,8);
for epoch = 1:851
    % use the transition matrix to prograte the state estimates
    x_minus = phi*x_est;

    % prograte the state estimates
    P_minus = phi*P_matrix*phi' + Q;

    % predict the ranges
    % calculate r_ej
    Sat_r_es_e = zeros(8,3);
    Sat_v_es_e = zeros(8,3);
    csv_Pseudo_ranges = csvread('Pseudo_ranges.csv');
    j = csv_Pseudo_ranges(1,2:size(csv_Pseudo_ranges,2));% j:satellite numbers
    % sat_r_es_e:ECEF satellite position, sat_v_es_e:ECEF satellite velocity
    for i = 1:size(j,2)
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(0.5*(epoch-1),j(i));
        Sat_r_es_e(i,:) = sat_r_es_e;%%r_ej
        Sat_v_es_e(i,:) = sat_v_es_e;%%v_ej
    end

    raj = zeros(8,1);
    for i = 1:8
        raj_tem = 0;
        for iteration = 1:10
            C = [1,omega_ie*raj_tem/c,0;-omega_ie*raj_tem/c,1,0;0,0,1];
            raj_tem = sqrt((C*Sat_r_es_e(i,:)' - x_minus(1:3))'*(C*Sat_r_es_e(i,:)' - x_minus(1:3)));
        end
        raj(i) = raj_tem;
    end

    % compute line-of-sight
    u = zeros(8,3);
    for i = 1:8
        C = [1,omega_ie*raj(i)/c,0;-omega_ie*raj(i)/c,1,0;0,0,1];
        u(i,:) = (C*Sat_r_es_e(i,:)' - x_minus(1:3))/raj(i);
    end

    % predict range rates
    der_raj = zeros(8,1);
    omega = [0 -omega_ie 0;omega_ie 0 0;0 0 0];
    for i = 1:8
        der_raj_tem = 0;
        C = [1,omega_ie*raj(i)/c,0;-omega_ie*raj(i)/c,1,0;0,0,1];
        for iteration = 1:10
            der_raj_tem = u(i,:)*(C*(Sat_v_es_e(i,:)' + omega*Sat_r_es_e(i,:)')-(x_minus(4:6)+omega*x_minus(1:3)));
        end
        der_raj(i) = der_raj_tem;
    end

    % compute the measurement matrix
    H_k = [-u,zeros(8,3),ones(8,1),zeros(8,1);
           zeros(8,3),-u,zeros(8,1),ones(8,1)];

    % compute the measurement noise convariance matrix
    sigma_rho = 10;
    sigma_r = 0.05;
    sigma_rho_matrix = sigma_rho^2*eye(8);
    sigma_r_matrix = sigma_r^2*eye(8);
    R_k = [sigma_rho_matrix,zeros(8,8);zeros(8,8),sigma_r_matrix];

    % compute the Kalman gain matrix
    K_k = P_minus*H_k'/(H_k*P_minus*H_k'+R_k);

    % formulate the measurement innovation vector
    % load measured pseudo-range
    csv_Pseudo_ranges = csvread('Pseudo_ranges.csv');
    rho = csv_Pseudo_ranges(epoch+1,2:size(csv_Pseudo_ranges,2))';
    % load measured pseudo-range rate
    csv_Pseudo_range_rates = csvread('Pseudo_range_rates.csv');
    der_rho = csv_Pseudo_range_rates(epoch+1,2:size(csv_Pseudo_range_rates,2))';
    delta_z = [rho-raj-x_minus(7)*ones(8,1);
               der_rho-der_raj-x_minus(8)*ones(8,1)];
    
    % detect outliers
    H_G = [-u,ones(8,1)]; 
    v = (H_G/(H_G'*H_G)*H_G'-eye(8))*delta_z(1:8);
    sigma_rho = 10;
    C_v = (eye(8) - H_G/(H_G'*H_G)*H_G')*sigma_rho^2;
    Threshold = 4;
    for number = 1:8
        tem = abs(v(number))- sqrt(C_v(number,number))*Threshold;
        if tem>0
        outlier_judge(epoch,number) = tem';
        end
    end
    if sum(outlier_judge(epoch,:) == 0) ~= 8
        [~,satellite] = max(outlier_judge(epoch,:));
%         satellite = 3;
        H_k(satellite+8,:)=[];
        H_k(satellite,:)=[];
        sigma_rho_matrix = sigma_rho^2*eye(7);
        sigma_r_matrix = sigma_r^2*eye(7);
        R_k = [sigma_rho_matrix,zeros(7,7);zeros(7,7),sigma_r_matrix];
        K_k = P_minus*H_k'/(H_k*P_minus*H_k'+R_k);
        delta_z(satellite+8,:)=[];
        delta_z(satellite,:)=[];
    end
    
    
    % update the state estimates
    x_pulse =  x_minus + K_k * delta_z;
    r_eb_e = x_pulse(1:3);
    v_eb_e = x_pulse(4:6);
    x_est = x_pulse;


    % error convariance matrix
    P_pulse = (eye(8) - K_k*H_k)*P_minus;
    P_matrix = P_pulse;

    % convert from ECEF to NED
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    latitude = L_b * rad_to_deg;
    longitude = lambda_b * rad_to_deg;
    position(epoch,:) = [latitude,longitude,h_b];
    velocity(epoch,:) = v_eb_n;
        
end