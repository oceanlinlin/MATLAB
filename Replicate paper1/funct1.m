clc
clear all
close all
%------f(x,y)=(sin(x)/x).*(sin(y)/y)-------
x = -10:0.02:10;
y = -10:0.02:10;
N = size(x,2);
for i=1:N
    for j=1:N
        z(i,j) = (sin(x(i))/x(i)).*(sin(y(j))/y(j));
    end
end
a = max(z);
mesh(x,y,z)
grid;
xlabel("x");
ylabel("y");
title("f(x,y)=(sin(x)/x)*(sin(y)/y)")


%函数在点(0,0)处可取得唯一最大极值1，同时周围分布着许多局部极值。