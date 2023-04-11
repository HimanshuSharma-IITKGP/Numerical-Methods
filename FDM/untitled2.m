close all

load = 0:200:1600;
r = 9.26/2000;
l = 0.120;
A = pi*r*r;
D1 = [0 1 2 3.5 5 7 9 10 12]*1e-3;
D2 = [0 1 3 5 6.5 8.5 10 11.5 13.5]*1e-3;
D_avg = (D1+D2)/2 ;

sigma = 1e-6*load/A;
strain = D_avg./l ;

plot(strain, sigma, '*')
coef = polyfit(strain, sigma, 1);
xfit = linspace(0, 150*1e-3, 1000);
yfit = polyval(coef, xfit);
hold on
plot(xfit, yfit, linewidth=2)

%
load2 = 1600:-200:0;
D1_2 = [12 8 7 5.5 3.5 2 1 0.5 0]*1e-3;
D2_2 = [13.5 10 9 7 5.5 4 2.5 1 0]*1e-3;
D_avg2 = (D1_2 + D2_2)/2 ;

sigma2 = 1e-6*load2/A;
strain2 = D_avg2./l ;

plot(strain2, sigma2, 'o')
coef2 = polyfit(strain2, sigma2, 1);
yfit2 = polyval(coef2, xfit);

plot(xfit, yfit2, linewidth=2)


