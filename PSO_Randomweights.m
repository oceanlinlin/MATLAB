%-------随机权重粒子群算法-----

format long;
%------初始化参数------
c1 = 2;
c2 = 2;
wmax = 0.9;
wmin = 0.4;
M = 100;%最大迭代次数
D = 2;%搜索空间维数
N = 50;%初始化种群数目
%函数式子没有给自变量的范围就不写
%速度没给具体的范围也先不写
Xmax = 10; %位置最大值
Xmin = -10; %位置最小值
Vmax = 10;  %速度最大值
Vmin = -10; %速度最小值

rande = 0.5; %随机权重方差

%------初始化种群的位置和速度------
for i = 1:N
    for j = 1:D
        x(i,j) = rand;
        v(i,j) = rand;
    end
end

%-----初始化个体最优位置与最优值-----
p = x;%第一代，最优位置即为初始位置
pbest = ones(N,1);%构建N*1矩阵，用于存储个体最优值
for i=1:N
    pbest(i) = fitness(x(i,:));  %各个粒子的最优适应度（最优值），因为是第一代，个体最优值即为函数本身的最优适应度
end

%------初始化全局最优位置与最优值------
gbest = inf;%初始设置该值为无穷大
g = zeros(1,D);%初始化最优位置矩阵，设置为全0矩阵
for i = 1:N
    if(pbest(i)<gbest)
        gbest = pbest(i);
        g = p(i,:);
    end
end
gb = ones(1,M);   %用于后续绘图，显示每次迭代的全局最优值

%----------------按照随机权重公式迭代直到满足精度或迭代次数--------------------
for t = 1:M
    %先设置随机权重法中的miu
    miu = wmin +(wmax-wmin)*rand;
    w = miu + rande * randn; %randn返回一个从标准正态分布中得到的随机标量
    for i=1:N
        %-------更新每个粒子所在的位置和速度------
        v(i,:) = w*v(i,:)+c1*rand*(p(i,:)-x(i,:))+c2*rand*(g-x(i,:));
        x(i,:) = x(i,:)+v(i,:);
        %-------更新个体最优位置和最优值-----
        if(fitness(x(i,:))<pbest(i))%若此时的适应函数值小于之前的最优值
            pbest(i) = fitness(x(i,:));
            p(i,:) = x(i,:);
        end
        %-----------更新全局最优位置和最优值-------------
        if(pbest(i)<gbest)
            gbest = pbest(i);
            g = p(i,:);
        end
        
    end
    %-----------记录历代全局最优值-------------
    gb(t) = gbest;
end

g  %显示全局最优的粒子位置
gb(end) %显示全局最优的粒子对应的最优值
figure
plot(gb)
xlabel('迭代次数');
ylabel("适应度值");
title("适应度进化曲线");
legend('RPSO');
 %-----------适应度函数-------------
function result = fitness(x)
result = x(1)^2+x(2)^2-x(1)*x(2)-10*x(1)-4*x(2)+60;
end