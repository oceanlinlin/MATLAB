%---- f(x,y)=(x+2y-7)^2+(2x+y-5)^2------

%绘制图像
clc
clear all
close all
[x,y] = meshgrid(-10:0.02:10,-10:0.02:10);
z = (x+2*y-7).^2+(2*x+y-5).^2;
figure;
mesh(x,y,z);
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
title('f(x,y) = (x+2y-7)^2+(2x+y-5)^2 ');


%标记出最大值点
maxVal = max(z(:));
[maxIndexX,maxIndexY] = find(z == maxVal);
for i = 1:length(maxIndexX)
    plot3(x(maxIndexX(i),maxIndexY(i)),y(maxIndexX(i),maxIndexY(i)), maxVal, 'r*','linewidth',2)
     text(x(maxIndexX(i),maxIndexY(i)),y(maxIndexX(i),maxIndexY(i)), maxVal, {['    X: ' num2str(x(maxIndexX(i),maxIndexY(i)))];['    Y: ' num2str(y(maxIndexX(i),maxIndexY(i)))];['    Z: ' num2str(maxVal)]})
    hold on
end

%标记出最小值点
minVal = min(z(:));
[minIndexX,minIndexY] = find(z == minVal);
for i = 1:length(minIndexX)
    plot3(x(minIndexX(i),minIndexY(i)),y(minIndexX(i),minIndexY(i)), minVal, 'r*','linewidth',2)
     text(x(minIndexX(i),minIndexY(i)),y(minIndexX(i),minIndexY(i)), minVal, {['    X: ' num2str(x(minIndexX(i),minIndexY(i)))];['    Y: ' num2str(y(minIndexX(i),minIndexY(i)))];['    Z: ' num2str(minVal)]})
    hold on
end