% Craft Temperature Cycle rev1.0 4/17/20
% Slugsat Science Experiment Subteam
% Tomohiro Shimada

%%
close all;
clear all;
clc;

t = linspace(0, 63072000, 63072000);
T = 5400;
f = 1/T;
D = 66;
temp_cycle = 62.5*square(2*pi*f.*t, D) + 22.5;

figure; plot(t, temp_cycle);
title("Simplified Max/Min Craft Temperature Cycle");
xlim([0 5400]); ylim([-100 120]);
xlabel("Time (seconds)"); ylabel("Craft Temperature [C]");
