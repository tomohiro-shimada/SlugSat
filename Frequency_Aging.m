% Frequency Aging Simulation rev1.0 4/10/20
% Slugsat Science Experiment Subteam
% Tomohiro Shimada

%%
close all;
clear all;
clc;

% time array
t = linspace(0,63072000,63072000);

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

% Plots
figure;
plot(t, a_t, 'g');
hold on
plot(t, b_t, 'b');
plot(t, c_t, 'y');
plot(t, Max, 'r--');
plot(t, Min, 'r--');
ylim([-2 2]);
xlim([1 63072000]);
title("Frequency Aging");
ylabel("Frequency Aging Deviation [ppm]");
xlabel("Time [seconds]");
legend("Stress Dominant", "Contamination Dominant", "Sum of Both", "Linear Case");