function results = Itegration(position,velocity,position_DR,velocity_DR,gyro_heading_deg)
%parameters
deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor

latitude_deg = position(:,1);
longitude_deg = position(:,2);
height = position(:,3);
v_N = velocity(:,1);
v_E = velocity(:,2);
latitude_rad = latitude_deg*deg_to_rad;
longitude_rad = longitude_deg*deg_to_rad;

DR_latitude_deg = position_DR(:,1);
DR_longitude_deg = position_DR(:,2);
DR_latitude_rad = DR_latitude_deg*deg_to_rad;
DR_longitude_rad = DR_longitude_deg*deg_to_rad;
DR_v_N_instant = velocity_DR(:,1);
DR_v_E_instant = velocity_DR(:,2);

%initialize P_0
x_pre = zeros(4,1);
%initialize P_0
sigma_v = 0.05;
sigma_r = 10;

times = size(position,1); %calculate the number of terms

[R_N,R_E]= Radii_of_curvature(latitude_rad(1));
P_pre = [sigma_v^2*eye(2),zeros(2,2);
         zeros(1,2),sigma_r^2/(R_N+height(1))^2,0;
         zeros(1,3),sigma_r^2/((R_E+height(1))^2*cos(latitude_rad(1))^2)];
 
results = zeros(351,6);
for i = 1:times-1

    %1.compute the transition matrix
    [R_N,R_E]= Radii_of_curvature(latitude_rad(i));
    tau = 0.5;
    phi_k_min_one = eye(4);
    phi_k_min_one(3,1) = tau/(R_N + height(i));
    phi_k_min_one(4,2) = tau/((R_E + height(i))*cos(latitude_rad(i)));

    %2.compute the system noise convariance matrix
    S_DR = 0.01;
    Q_k_min_one = zeros(4,4);
    Q_k_min_one(1,1) = S_DR*tau;
    Q_k_min_one(2,2) = Q_k_min_one(1,1);
    Q_k_min_one(3,3) = S_DR*tau^3/(R_N+height(i))^2/3;
    Q_k_min_one(4,4) = S_DR*tau^3/((R_E+height(i))^2*cos(latitude_rad(i))^2)/3;
    Q_k_min_one(1,3) = S_DR*tau^2/(R_N+height(i))/2;
    Q_k_min_one(3,1) = Q_k_min_one(1,3);
    Q_k_min_one(2,4) = S_DR*tau^2/((R_E+height(i))*cos(latitude_rad(i)))/2;
    Q_k_min_one(4,2) = Q_k_min_one(2,4);

    %3.compute state estimate
    x_k_minus = phi_k_min_one*x_pre;

    %4.compute the error convariance matrix
    P_k_minus = phi_k_min_one*P_pre*phi_k_min_one' + Q_k_min_one;

    %5.measurement matrix
    H_k = [zeros(2,2),-eye(2);-eye(2),zeros(2,2)];

    %6.measurement noise convariance matrix
    sigma_Gr = 5.2;
    sigma_Gv = 0.02;
    R_k = zeros(4,4);
    [R_N,R_E]= Radii_of_curvature(latitude_rad(i+1)); %new R_N and R_E
    R_k(1,1) = sigma_Gr^2/(R_N + height(i+1))^2;
    R_k(2,2) = sigma_Gr^2/((R_E + height(i+1))^2*cos(latitude_rad(i+1))^2);
    R_k(3,3) = sigma_Gv^2;
    R_k(4,4) = sigma_Gv^2;

    %7.Kalman gain matrix
    K_k = P_k_minus*H_k'/(H_k*P_k_minus*H_k'+R_k);

    %8.measurment innovation vector
    delta_z_k_minus = [latitude_rad(i+1)-DR_latitude_rad(i+1);
                       longitude_rad(i+1)-DR_longitude_rad(i+1);
                       v_N(i+1) - DR_v_N_instant(i+1);
                       v_E(i+1) - DR_v_E_instant(i+1)] - H_k*x_k_minus;

    %9.update the state estimates
    x_k_plus = x_k_minus + K_k*delta_z_k_minus;
    x_pre = x_k_plus;
    
    %10.update the error convariance matrix
    P_k_plus = (eye(4) - K_k*H_k)*P_k_minus;
    P_pre = P_k_plus;

    %11.results
    v_N_corrected = DR_v_N_instant(i+1) - x_k_plus(1);
    v_E_corrected = DR_v_E_instant(i+1) - x_k_plus(2);
    
    latitude_corrected_rad = DR_latitude_rad(i+1) -  x_k_plus(3);
    longitude_corrected_rad = DR_longitude_rad(i+1) -  x_k_plus(4);
    latitude_corrected_deg = latitude_corrected_rad * rad_to_deg;
    longitude_corrected_deg = longitude_corrected_rad * rad_to_deg;

    results(i+1,1) = 0.5*i;
    results(i+1,2) = latitude_corrected_deg;
    results(i+1,3) = longitude_corrected_deg;
    results(i+1,4) = v_N_corrected;
    results(i+1,5) = v_E_corrected;
end

results(1,1) = 0;
results(1,2) = latitude_deg(1);
results(1,3) = longitude_deg(1);
results(1,4) = v_N(1);
results(1,5) = v_E(1);
results(:,6) = gyro_heading_deg;

%plot position figure
plot(results(:,3),results(:,2));hold on;
plot(longitude_deg,latitude_deg);hold on;
plot(DR_longitude_deg,DR_latitude_deg);hold on;
xlabel('Longitude in degree');
ylabel('Latitude in degree');
legend('Integration','GNSS','DR')
figure;

%plot velocity figure
plot(results(:,1),results(:,4));hold on;
plot(results(:,1),results(:,5));hold on;
xlabel('Time');
ylabel('Velocity in meters per second');
legend('North velocity','East velocity')
figure;

%plot heading
plot(results(:,1),results(:,6));hold on;
xlabel('Time');
ylabel('Heading in degrees');

