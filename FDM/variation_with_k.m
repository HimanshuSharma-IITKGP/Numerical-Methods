clear;
clc;
close all;

for k = 1:2:9
    figure;

    C = diffusion(1.2, 0, 1, 11, 1e-6, 1e-12, k, 0.004, 1);
    disp(C(1))
    ylim([0, 1.8]) ;
    pause(0.1) ;
end