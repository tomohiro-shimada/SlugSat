
close all;
clear all;
clc;
% time array
t = linspace(0,63072000,63072000);
temp = randi([-40 85], 1,length(t))+zeros(1, length(t));
[stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp);

% Frequency Aging Variables/Scenarios
[a_t, b_t, c_t, Max, Min] = frequency_aging(t);

Temp1 = a_t + stab_tol_values_temp;
Temp2 = b_t + stab_tol_values_temp;
Temp3 = c_t + stab_tol_values_temp;
Temp4 = Min + stab_tol_values_temp;
Temp5 = Max + stab_tol_values_temp;
total_drift_sec_a = trapz(Temp1)/(10^6)
total_drift_sec_b = trapz(Temp2)/(10^6)
total_drift_sec_c = trapz(Temp3)/(10^6)
total_drift_sec_Min = trapz(Temp4)/(10^6)
total_drift_sec_Max = trapz(Temp5)/(10^6)

% Plots
figure;
plot(t, a_t, 'r--');
hold on
plot(t, Temp1, 'b');
ylim([-1 1]);
xlim([1 50]);
title("Frequency Aging/Tolerance/Stability");
ylabel("Frequency Aging Deviation [ppm]");
xlabel("Time [seconds]");
legend("Aging", "Total");

function [a_t, b_t, c_t, Max, Min] = frequency_aging(t)
% Frequency Aging with Stress as Dominating Factor
x_s = 3.017*10^-2;
a_t = x_s*log(0.5*t + 1);

% Frequency Aging with Contamination as Dominating Factor
x_c = 4.115*10^-2;
b_t = -x_c*log(0.006*t +1);

% The Sum of Both Equations
c_t = a_t + b_t;

% Linear Case 
Max = 1/length(t).* t;
Min = -1/length(t).* t;
end

function [stab_tol_values_rand, stab_tol_values_temp] = stab_tol(temp)

T = 120;
f = 1/T;
ii = 1;
stab_tol_values_temp = zeros(1, length(temp));
while ii <= length(temp)
    stab_tol_values_temp(ii) = 0.28*sin(2*pi*f.*(temp(ii) + 35)) + 0.005*randi([-100 100]);
    ii = ii + 1;
end

stab_tol_values_rand = 0.0078.*randi([-100 100], 1, length(temp));
end