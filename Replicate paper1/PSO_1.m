%------标准粒子群算法求函数极值------
%------给定初始化条件-----------
clear all;
close all;
clc;
N = 20; %群体粒子个数
D = 2;  %粒子搜索空间维数
T = 500; %最大迭代次数
c1 = 1.5; %学习因子1
c2 = 1.7; %学习因子2
w = 1.2;  %惯性权重
Xmax = 10; %位置最大值
Xmin = -10; %位置最小值
Vmax = 100;  %速度最大值
Vmin = -100; %速度最小值

format long;

%------初始化种群的个体（可以在此处限定位置和速度的范围）------------
x = rand(N,D) * (Xmax-Xmin) + Xmin;
v = rand(N,D) * (Vmax-Vmin) + Vmin;


%------初始化个体最优位置和最优值----------------------

p = x;  %因为是第一代，所以初始时个体最优位置p即为初始位置
pbest = ones(N,1);  %构建N*1矩阵，用于存储个体最优值
for i=1:N
    pbest(i) = func1(x(i,:));  %各个粒子的最优适应度（最优值），因为是第一代，个体最优值即为函数本身的最优适应度
end

%------初始化全局最优位置和全局最优值------------
g = ones(1,D); %构建1*N矩阵，由于存储全局最优位置
gbest = 0;  %初始时，全局最优值设为0，之后比较时选取大于号用于求取函数最大值
%若初始时候gbest取值为inf,之后比较时用小于号来求取最小值
for i = 1:N
    if(pbest(i)>gbest)  %大于符号用于求取最大值，小于符号用于求取最小值
        gbest = pbest(i);  %此时的全局最优值为第i个粒子的最优值
        g = p(i,:);    %此时的全局最优位置为粒子矩阵p中第i个粒子对应的初始位置
    end
end
gb = ones(1,T);   %用于后续绘图，显示每次迭代的全局最优值


%----------------按照公式迭代直到满足精度或迭代次数----------------------
for i=1:T
    for j=1:N
        %-----------更新个体最优位置和最优值-------------
        if(func1(x(j,:))>pbest(j)) %若此时的适应函数值大于之前的最优值，目的是求最大值；若为小于符号，求的是最小值
            pbest(j) = func1(x(j,:));
            p(j,:) = x(j,:);   %此时的个体最优位置等于取当前最优值时候的粒子位置
        end
        %-----------更新全局最优位置和最优值-------------
        if(pbest(j)>gbest) %若为大于符号求最大值，小于符号求最小值
            gbest = pbest(j); %此时的全局最优值等于当前粒子的个体最优值
            g = p(j,:);
        end
        %-----------更新位置和速度值-------------
        v(j,:) = w * v(j,:) + c1 * rand *(p(j,:)-x(j,:)) + c2 * rand *(g - x(j,:));
        x(j,:) = x(j,:) + v(j,:);
        %-----------边界条件处理-------------
        for ii = 1:D
            if(v(j,ii) > Vmax || v(j,ii) < Vmin)
                v(j,ii) = rand * (Vmax-Vmin) + Vmin;
            end
            if(x(j,ii) > Xmax || x(j,ii) < Xmin)
                x(j,ii) = rand *(Xmax - Xmin) + Xmin;
            end
        end
    end
    %-----------记录历代全局最优值-------------
    gb(i) = gbest;
end
g  %显示全局最优的粒子位置
gb(end) %显示全局最优的粒子对应的最优值
hf2 = figure;
%figure;
plot(gb)
xlabel('迭代次数');
ylabel("适应度值");
title("适应度进化曲线");
legend('PSO');
savefig(hf2,"PSO曲线.fig");
 %-----------适应度函数-------------
function result = func1(x)
result = (sin(x(1))./x(1)).*(sin(x(2))./x(2));
end

