%-----f(x,y) = 20*exp(-0.2*sqrt(0.5*(x^2+y^2)))
%+exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))-exp(1)-20-----

clc
clear all
close all
x = -30:0.02:30;
y = -30:0.02:30;
N = size(x,2);
for i=1:N
    for j=1:N
        z(i,j) = 20*exp(-0.2*sqrt((x(i)^2+y(j)^2)/2)) + exp((cos(2*pi*x(i))+cos(2*pi*y(j)))/2) - exp(1) - 20;
    end
end
a = max(z);
mesh(x,y,z)
grid
xlabel("x");
ylabel("y");
title("f(x,y)=20*exp(-0.2*sqrt(0.5*(x^2+y^2)))+exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))-exp(1)-20")