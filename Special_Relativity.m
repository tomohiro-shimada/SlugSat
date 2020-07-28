% Special Relativity Simulations rev1.0 4/24/20
% Slugsat Science Experiment Subteam
% Tomohiro Shimada

%%
close all;
clear all;
clc; 

% Variable Definition
t = linspace(0,63072000,63072000);
r_earth = 6378*10^3;
r_sat = 20000*10^3;
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


figure;
plot(t, Norm_SR_TD);
hold on;
plot(t, Norm_GR_TD);
xlim([0 63072000]);
title("Relativistic Effects on GPS Clock");
xlabel("Time (seconds)");
ylabel("Drift (Seconds)");
legend("Special Relativity", "General Relativity");

