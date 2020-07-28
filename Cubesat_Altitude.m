% Simulated Craft Altitiude rev1.0 4/17/20
% Slugsat Science Experiment Subteam
% Tomohiro Shimada

%%
close all;
clear all;
clc; 

altitudes = [100 200 300 400 500 600 700 800 900];
T = 5400;
t = 63072000;
craft_altitude_period_n = randi([1 9], 1, t/T).*100;
period_time = linspace(0,t/T, t/T);

lifetime_altitude = zeros(1, t);
ii = 1;
while ii < t/T
    lifetime_altitude(ii:ii+5399) = (craft_altitude_period_n(ii+1)-craft_altitude_period_n(ii))/(period_time(ii+1)-period_time(ii))+craft_altitude_period_n(ii);
    ii = ii + 1;
end

figure;
plot(t, lifetime_altitude);
ylim([50 1000]);


        
        
        
