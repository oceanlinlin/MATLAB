%-----f(x,y) = 20+(x^2-10*cos(2*pi*x))+(y^2-10*cos(2*pi*y))-----

clc
clear all
close all
x = -512:1:512;
y = -512:1:512;
N = size(x,2);
for i=1:N
    for j=1:N
        z(i,j) = 20+(x(i)^2-10*cos(2*pi*x(i)))+(y(j)^2-10*cos(2*pi*y(j)));
    end
end
a = min(z);
mesh(x,y,z)
grid
xlabel("x");
ylabel("y");
title("f(x,y)=20+(x^2-10*cos(2*pi*x))+(y^2-10*cos(2*pi*y))");