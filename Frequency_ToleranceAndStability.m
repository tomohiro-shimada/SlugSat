% Frequency Tolerance and Stability Simulation rev1.0 4/10/20
% Slugsat Science Experiment Subteam
% Tomohiro Shimada

%%
close all;
clear all;
clc;

% Operating Temperature Range 
temp = linspace(-40, 85, 126);

% Frequency Stability Curve
T = 120;
f = 1/T;
freq_stab_curve = 0.28*sin(2*pi*f.*(temp + 35));

% Frequency Stability + Frequency Tolerance
freq_drift_max = freq_stab_curve + 0.5;
freq_drift_min = freq_stab_curve - 0.5;

% Actual Total Short term Stability
freq_drift_temp = zeros(1, length(temp));
ii = 1;
while ii <= length(temp)
    freq_drift_temp(ii) = freq_stab_curve(ii) + 0.01*randi([-50 50]); 
    ii = ii +1;
end

figure;
plot(temp, freq_stab_curve, 'b--');
hold on
plot(temp, freq_drift_max, 'r--');
plot(temp, freq_drift_temp, 'g');
plot(temp, freq_drift_min, 'r--');
ylim([-1 1]);
xlim([-40 85]);
title("Frequency Stability + Tolerance, with Temperature Dependance");
ylabel('Short Term Drift [ppm]')
xlabel('Temperature [C]');
legend('Frequency Stabillity', 'Total Frequency Drift Min/Max', 'Total Short Term Stability');

total_min = 0.78*(ones(1,length(temp)));
total_max = -0.78*(ones(1,length(temp)));

freq_drift_rand = zeros(1, length(temp));
ii = 1;
while ii <= length(temp)
    freq_drift_rand(ii) = 0.0078*randi([-100 100]);
    ii = ii+1;
end

figure;
plot(temp, freq_drift_rand, 'g');
hold on
plot(temp, total_max, 'r--');
plot(temp, total_min, 'r--');
ylim([-1 1]);
xlim([-40 85]);
title("Frequency Stability + Tolerance, with No Temperature Dependance");
ylabel('Short Term Drift [ppm]')
xlabel('Temperature [C]');
legend('Total Short Term Stability', 'Total Frequency Drift Min/Max');
