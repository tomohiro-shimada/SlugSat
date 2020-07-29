% Tomohiro Shimada
% Slugsat Science Experiment SubTeam

%Main File for RTC Simulations

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

% time array
t = linspace(0,63072000,63072000);

% Craft Temperature Cycle
T = 5400;       %Period of Temperature Change
f = 1/T;
D = 66;     %Duty Cycle of Temperature Square Wave
temp = 62.5*square(2*pi*f.*t, D) + 22.5;  %Craft Temperature over 2 years

% Special/General Relativity Variables
r_earth = 6378*10^3; %Radius of Earth
r_sat = input('Enter operating altitude of Cubesat in km: ') * 10^3;    %Operating Altitude of Cubesat

% Time Dilation Function Call
total_TD = TD_calculation(r_earth, r_sat, t);

% Frequency Aging Variables/Scenarios Function Call
f_a = input('Enter frequency aging value of crystal oscillator in ppm/year: ');
[a_t, b_t, c_t, Max, Min] = frequency_aging(t, f_a);

% Frequency Stability/Tolerance  function Call
stab = input('Enter frequency stability value of crystal oscillator in ppm: ');
tol  = input('Enter frequency tolerance value of crystal oscillator in ppm: ');
[stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp, stab, tol);

% Frequency characteristics of crystal translated to drift in seconds
Temp1 = a_t + stab_tol_values_temp;
Temp2 = b_t + stab_tol_values_temp;
Temp3 = c_t + stab_tol_values_temp;
Temp4 = Min + stab_tol_values_temp;
Temp5 = Max + stab_tol_values_temp;

% Integration of drift as a function of time (i.e. gets total time at the
% end)
ii = 2;
running_total_1 = zeros(1, length(t));
running_total_2 = zeros(1, length(t));
running_total_3 = zeros(1, length(t));
running_total_4 = zeros(1, length(t));
running_total_5 = zeros(1, length(t));
while ii <= length(t)
   running_total_1(ii) = running_total_1(ii-1) + Temp1(ii);
   running_total_2(ii) = running_total_2(ii-1) + Temp2(ii);
   running_total_3(ii) = running_total_3(ii-1) + Temp3(ii);
   running_total_4(ii) = running_total_4(ii-1) + Temp4(ii);
   running_total_5(ii) = running_total_5(ii-1) + Temp5(ii);
   ii = ii+1;
end
running_total_1 = running_total_1./10^6 + total_TD;
running_total_2 = running_total_2./10^6 + total_TD;
running_total_3 = running_total_3./10^6 + total_TD;
running_total_4 = running_total_4./10^6 + total_TD;
running_total_5 = running_total_5./10^6 + total_TD;

%%
% Simulated downlink (1 per day)
t_sampled = t(1:86400:end);
downlinked_data_1 = running_total_1(1:86400:end);
downlinked_data_2 = running_total_2(1:86400:end);
downlinked_data_3 = running_total_3(1:86400:end);
downlinked_data_4 = running_total_4(1:86400:end);
downlinked_data_5 = running_total_5(1:86400:end);

% Adding propagation delay error to downlink times
ii = 1;
while ii <= length(downlinked_data_1)
   downlinked_data_1(ii) = downlinked_data_1(ii) + 5.3*10^-6*randi([-10 10]);
   downlinked_data_2(ii) = downlinked_data_2(ii) + 5.3*10^-6*randi([-10 10]);
   downlinked_data_3(ii) = downlinked_data_3(ii) + 5.3*10^-6*randi([-10 10]);
   downlinked_data_4(ii) = downlinked_data_4(ii) + 5.3*10^-6*randi([-10 10]);
   downlinked_data_5(ii) = downlinked_data_5(ii) + 5.3*10^-6*randi([-10 10]);
   ii = ii + 1;
end

interpolated_time_1 = interp1(t_sampled, downlinked_data_1, t);
interpolated_time_2 = interp1(t_sampled, downlinked_data_2, t);
interpolated_time_3 = interp1(t_sampled, downlinked_data_3, t);
interpolated_time_4 = interp1(t_sampled, downlinked_data_4, t);
interpolated_time_5 = interp1(t_sampled, downlinked_data_5, t);


figure;
plot(t_sampled, downlinked_data_2, 'ko');
hold on;
plot(t, running_total_2, 'b');
xlim([0 86400*10]);
ylabel("Drift (seconds)");
xlabel("Time During Orbit (seconds)");
title("Simulated Downlinked Points of Board Time Drift")
legend("Downlink", "On Board Clock Drift");

% actual time - interpolated time
offset_1 = running_total_1 - interpolated_time_1;
max_offset_1 = max(abs(offset_1));
offset_2 = running_total_2 - interpolated_time_2;
max_offset_2 = max(abs(offset_2));
offset_3 = running_total_3 - interpolated_time_3;
max_offset_3 = max(abs(offset_3));
offset_4 = running_total_4 - interpolated_time_4;
max_offset_4 = max(abs(offset_4));
offset_5 = running_total_5 - interpolated_time_5;
max_offset_5 = max(abs(offset_5));

%% 
% Statistical Analysis
% Finding Equations of Regression Lines and Correlation Coefficients
y_intercept = zeros(5, 730);
y_intercept(1,:) = downlinked_data_1;
y_intercept(2,:) = downlinked_data_2;
y_intercept(3,:) = downlinked_data_3;
y_intercept(4,:) = downlinked_data_4;
y_intercept(5,:) = downlinked_data_5;

slope = zeros(5, 729);
ii = 1;
while ii<=729
    slope(1, ii) = (y_intercept(1,ii+1) - y_intercept(1,ii))/86400;
    slope(2, ii) = (y_intercept(2,ii+1) - y_intercept(2,ii))/86400;
    slope(3, ii) = (y_intercept(3,ii+1) - y_intercept(3,ii))/86400;
    slope(4, ii) = (y_intercept(4,ii+1) - y_intercept(4,ii))/86400;
    slope(5, ii) = (y_intercept(5,ii+1) - y_intercept(5,ii))/86400;
    ii = ii + 1;
end
xlswrite('slope_of_regression_lines.xlsx', slope);
y_intercept(:, 730) = [];

xlswrite('y_intercept_of_regression_lines.xlsx', y_intercept);

r1 = [];
r2 = [];
r3 = [];
r4 = [];
r5 = [];
ii = 1;
while ii <= 729
    n1 = (ii-1)*86400 + 1;
    n2 = ii*86400;
    temp1 = corrcoef(running_total_1(1, n1:n2),interpolated_time_1(1,n1:n2));
    r_temp1 = temp1(2,1); 
    r1 = [r1 r_temp1];
    
    temp2 = corrcoef(running_total_2(1, n1:n2),interpolated_time_2(1,n1:n2));
    r_temp2 = temp2(2,1); 
    r2 = [r2 r_temp2];
 
    temp3 = corrcoef(running_total_3(1, n1:n2),interpolated_time_3(1,n1:n2));
    r_temp3 = temp3(2,1); 
    r3 = [r3 r_temp3];
    
    temp4 = corrcoef(running_total_4(1, n1:n2),interpolated_time_4(1,n1:n2));
    r_temp4 = temp4(2,1); 
    r4 = [r4 r_temp4];
    
    temp5 = corrcoef(running_total_5(1, n1:n2),interpolated_time_5(1,n1:n2));
    r_temp5 = temp5(2,1); 
    r5 = [r5 r_temp5];
    
    ii = ii + 1;
end
%%
% Plots
figure;
plot(t, running_total_1, 'b');
hold on
plot(t, interpolated_time_1, 'r--');
plot(t_sampled, downlinked_data_1, 'ko');
% plot(t, running_total_2, 'g');
% plot(t, interpolated_time_2, 'g--');
% plot(t, running_total_3, 'k');
% plot(t, interpolated_time_3, 'k--');
% plot(t, running_total_4, 'r');
% plot(t, interpolated_time_4, 'r--');
% plot(t, running_total_5, 'r');
% plot(t, interpolated_time_5, 'r--');

%ylim([-0.003 0.003]);
xlim([0 86401*2]);
title("Actual and Approximated On Board Clock Drift vs Time");
ylabel("Drfit [Seconds]");
xlabel("Time [seconds]");
legend("Actual Drift", "Approximated Drift", "Downlink");

%%
% Functions

% Finding Frequency Aging Scenarios
function [a_t, b_t, c_t, Max, Min] = frequency_aging(t, f_a)

% Frequency Aging with Stress as Dominating Factor
x_s = f_a/(log(0.5*31536000+1));
a_t = -x_s*log(0.5*t + 1);

% Frequency Aging with Contamination as Dominating Factor
x_c = f_a/(log(0.0006*31536000 + 1));
b_t = x_c*log(0.006*t +1);

% The Sum of Both Equations
c_t = a_t + b_t;

% Linear Min/Max Cases
Max = -(f_a/31536000).* t;
Min = (f_a/31536000).* t;
end

% Finding Stability and Tolerance 
function [stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp,stab,tol)
% Variables
T = 120;
f = 1/T;

tol1 = tol/100;

%Temperature dependent Frequency stability
ii = 1;
stab_tol_values_temp = zeros(1, length(temp));
while ii <= length(temp)
    stab_tol_values_temp(ii) = -stab*sin(2*pi*f.*(temp(ii) + 35)) + tol1*randi([-100 100]);
    ii = ii + 1;
end

% Assuming Stability is not dependent on temperature
stab_tol_values_rand = 0.0078.*randi([-100 100], 1, length(temp));
end

% Time Dilation Function
function total_TD = TD_calculation(r_earth, r_sat, t)
G = 6.67*10^-11;
M = 5.98*10^24;
c = 3.0*10^8;

%Special Relativity
v = sqrt(G*M/(r_earth + r_sat)); % Finding Velocity 
TDF = 1/(sqrt(1-(v^2/c^2)));     % Time Dilation Factor
SR_TD = t*TDF;               % Time dilation comparason b/w cubesat and GS
Norm_SR_TD = SR_TD - t;      % resulting time difference

%General Relativity
GR_TD = (1-G*M/c^2*(1/r_earth-1/(r_earth+r_sat)))*t; % Time dilation comparason b/w cubesat and GS
Norm_GR_TD = GR_TD - t;        % Resulting Time Difference

total_TD = Norm_SR_TD - Norm_GR_TD;   % Total Time Difference between both effects
end