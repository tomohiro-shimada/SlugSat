close all;
clear all;
clc;

T = 86400;
t = linspace(0, 63072000, 63072000);
y = 10*sin(2*pi/(T)*t);

t_sampled = t(1:T:end);
y_sampled = y(1:T:end);

y_interp = interp1(t_sampled,y_sampled,t);

%%
figure;
plot(t, y, 'r-');
hold on;
plot(t, y_interp, 'b');
xlim([0 86400]);
