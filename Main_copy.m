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
r_earth = 6378*10^3;
r_sat = 600*10^3;

% Time Dilation
total_TD = TD_calculation(r_earth, r_sat, t);

% Frequency Stability/Tolerance 
[stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp);

% Frequency Aging Variables/Scenarios
[a_t, b_t, c_t, Max, Min] = frequency_aging(t);

Temp1 = a_t + stab_tol_values_temp;
Temp2 = b_t + stab_tol_values_temp;
Temp3 = c_t + stab_tol_values_temp;
Temp4 = Min + stab_tol_values_temp;
Temp5 = Max + stab_tol_values_temp;

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
interpolated_time_1 = interp1(t_sampled, downlinked_data_1, t);
downlinked_data_2 = running_total_2(1:86400:end);
interpolated_time_2 = interp1(t_sampled, downlinked_data_2, t);
downlinked_data_3 = running_total_3(1:86400:end);
interpolated_time_3 = interp1(t_sampled, downlinked_data_3, t);
downlinked_data_4 = running_total_4(1:86400:end);
interpolated_time_4 = interp1(t_sampled, downlinked_data_4, t);
downlinked_data_5 = running_total_5(1:86400:end);
interpolated_time_5 = interp1(t_sampled, downlinked_data_5, t);

figure;
plot(t_sampled, downlinked_data_2, 'ro');
xlim([0 86400*10]);
ylabel("Drift (seconds)");
xlabel("time in orbit (seconds)");
title("Simulated Downlinked On Board Time Drift")
legend("Downlinked Time");

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
% Plots
figure;
plot(t, running_total_1, 'b');
hold on
plot(t, interpolated_time_1, 'b--');
plot(t, running_total_2, 'g');
plot(t, interpolated_time_2, 'g--');
plot(t, running_total_3, 'k');
plot(t, interpolated_time_3, 'k--');
plot(t, running_total_4, 'r');
plot(t, interpolated_time_4, 'r--');
plot(t, running_total_5, 'r');
plot(t, interpolated_time_5, 'r--');

%ylim([-0.003 0.003]);
xlim([0 86400*2]);
title("Frequency Aging/Tolerance/Stability");
ylabel("Drfit in Seconds");
xlabel("Time [seconds]");


%%
% Functions

function [a_t, b_t, c_t, Max, Min] = frequency_aging(t)
% Frequency Aging with Stress as Dominating Factor
x_s = 3.017*10^-2;
a_t = -x_s*log(0.5*t + 1);

% Frequency Aging with Contamination as Dominating Factor
x_c = 4.115*10^-2;
b_t = x_c*log(0.006*t +1);

% The Sum of Both Equations
c_t = a_t + b_t;

% Linear Case 
Max = -1/length(t).* t;
Min = 1/length(t).* t;
end

function [stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp)

T = 120;
f = 1/T;
ii = 1;
stab_tol_values_temp = zeros(1, length(temp));
while ii <= length(temp)
    stab_tol_values_temp(ii) = -0.28*sin(2*pi*f.*(temp(ii) + 35)) + 0.015*randi([-100 100]);
    ii = ii + 1;
end

stab_tol_values_rand = 0.0078.*randi([-100 100], 1, length(temp));
end

function total_TD = TD_calculation(r_earth, r_sat, t)
G = 6.67*10^-11;
M = 5.98*10^24;
c = 3.0*10^8;

%Special Relativity
v = sqrt(G*M/(r_earth + r_sat));
TDF = 1/(sqrt(1-(v^2/c^2)));
SR_TD = t*TDF;
Norm_SR_TD = SR_TD - t;

%General Relativity
GR_TD = (1-G*M/c^2*(1/r_earth-1/(r_earth+r_sat)))*t;
Norm_GR_TD = GR_TD - t;

total_TD = Norm_SR_TD - Norm_GR_TD;
end