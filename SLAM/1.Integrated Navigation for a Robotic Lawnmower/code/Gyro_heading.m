function gyro_heading_deg = Gyro_heading()
% parameters
deg_to_rad = 0.01745329252; 
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor

Dead_reckoning = csvread('Dead_reckoning.csv');
gyro_angular_rate = Dead_reckoning(:,6);
magnetic_heading_deg = Dead_reckoning(:,7);
magnetic_heading_rad = magnetic_heading_deg*deg_to_rad;

% calculate the gyroscope headling
gyro_heading_deg = zeros(851,1);
gyro_heading_deg(1) = magnetic_heading_deg(1);
for i = 2:851
    gyro_heading_deg(i) = gyro_heading_deg(i-1) + (gyro_angular_rate(i)*0.5*rad_to_deg);
end
gyro_heading_rad = gyro_heading_deg * deg_to_rad;

% update heading error
x_pre = [0;0];
phi_k_minus_1 = [1,0.5;0,1];
P_k_minus_1 = [0.01^2,0;0,(1*deg_to_rad)^2];
%P_k_minus_1 = [1,0;1,0];
S_rg = 3*10^-6;
S_bgd = (1*deg_to_rad)^2;
%S_bgd = 1*10^-12;
H_k = [-1,0];

% calculate heading in each epoch
for t = 2:851
    Q_k_minus_1 = [S_rg*0.5+S_bgd*0.5^3/3,S_bgd*0.5^2/2;S_bgd*0.5^2/2,S_bgd*0.5];
    x_k_minus = phi_k_minus_1 * x_pre;
    P_k_minus = phi_k_minus_1 * P_k_minus_1 * phi_k_minus_1'+ Q_k_minus_1;
    delta_z_k_minus = (magnetic_heading_rad(t) - gyro_heading_rad(t)) - H_k*x_k_minus;
    R_k = (4*deg_to_rad)^2;
    K_k = P_k_minus*H_k'/(H_k * P_k_minus * H_k' + R_k);
    % update x
    x_pre = x_k_minus + K_k * delta_z_k_minus;
    % update P
    P_k_minus_1 = (eye(2) - K_k*H_k) * P_k_minus;
    gyro_heading_deg(t) = (gyro_heading_rad(t) - x_pre(1))*rad_to_deg;
end

