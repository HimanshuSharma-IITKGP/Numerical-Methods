
clear;
close all ;
clc; 

for D = 0.1:0.1:0.5
    figure ;
    C = diffusion(D, 0, 1, 51, 1e-6, 1e-12, 1, 0.0004, 0.1);
    pause(0.3) ;
end