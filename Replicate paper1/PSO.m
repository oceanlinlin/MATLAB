function [PSO_pbestloc,PSO_gbestloc,PSO_Best,PSO_end] = PSO(N,D,T,Xmax,Xmin,Vmax,Vmin)
c1 = 1.5;%学习因子
c2 = 1.7;
w = 1.2;%权重系数

%----初始化种群的速度和位置----
x = rand(N,D) * (Xmax-Xmin) + Xmin;
v = rand(N,D) * (Vmax-Vmin) + Vmin;

%------初始化个体最优位置和最优值----------------------

PSO_pbestloc = x;  %因为是第一代，所以初始时个体最优位置p即为初始位置
PSO_pbest = ones(N,1);  %构建N*1矩阵，用于存储个体最优值
for i=1:N
    PSO_pbest(i) = func1(x(i,:));  %各个粒子的最优适应度（最优值），因为是第一代，个体最优值即为函数本身的最优适应度
end

%------初始化全局最优位置和全局最优值------------
PSO_gbestloc = ones(1,D); %构建1*N矩阵，由于存储全局最优位置
PSO_Best = inf;  %初始时，全局最优值设为0，之后比较时选取大于号用于求取函数最大值
%论文中f1函数求最大值时候可以设置PSO_Best为0;f2函数得设置其值初始为负无穷，或者更直接一点，求最大值时候全设置为-inf
%论文中f3函数有最小值，可设置初始全局最优值为inf,下面对应的＞要改为＜
%若初始时候gbest取值为inf,之后比较时用小于号来求取最小值
for i = 1:N
    if(PSO_pbest(i)<PSO_Best)  %大于符号用于求取最大值，小于符号用于求取最小值
        PSO_Best = PSO_pbest(i);  %此时的全局最优值为第i个粒子的最优值
        PSO_gbestloc = PSO_pbest(i,:);    %此时的全局最优位置为粒子矩阵p中第i个粒子对应的初始位置
    end
end
gb = ones(1,T);   %用于后续绘图，显示每次迭代的全局最优值



%----------------按照公式迭代直到满足精度或迭代次数----------------------
for i=1:T
    for j=1:N
        %-----------更新个体最优位置和最优值-------------
        if(func1(x(j,:))<PSO_pbest(j)) %若此时的适应函数值大于之前的最优值，目的是求最大值；若为小于符号，求的是最小值
            PSO_pbest(j) = func1(x(j,:));
            PSO_pbestloc(j,:) = x(j,:);   %此时的个体最优位置等于取当前最优值时候的粒子位置
        end
        %-----------更新全局最优位置和最优值-------------
        if(PSO_pbest(j)<PSO_Best) %若为大于符号求最大值，小于符号求最小值
            PSO_Best = PSO_pbest(j); %此时的全局最优值等于当前粒子的个体最优值
            PSO_gbestloc = PSO_pbestloc(j,:);
        end
        %-----------更新位置和速度值-------------
        v(j,:) = w * v(j,:) + c1 * rand *(PSO_pbestloc(j,:)-x(j,:)) + c2 * rand *(PSO_gbestloc - x(j,:));
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
    gb(i) = PSO_Best;
end
PSO_gbestloc;  %显示全局最优的粒子位置
PSO_Best = gb; %显示全局最优的粒子对应的最优值
PSO_end = gb(end);
PSO_pbestloc;
%hf2 = figure;
%figure;
%plot(gb,'--k','LineWidth',2)
%xlabel('迭代次数');
%ylabel("适应度值");
%title("适应度进化曲线");
%legend('PSO');
%savefig(hf2,"PSO曲线.fig");


end

function result  = func1(x)
result = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
%result = x(1)^2 + x(2)^2;
%result = 20+(x(1)^2-10*cos(2*pi*x(1)))+(x(2)^2-10*cos(2*pi*x(2)));
%result = 20*exp(-0.2*sqrt((x(1)^2+x(2)^2)/2)) + exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2) - exp(1) - 20;
%result = (sin(x(1))/x(1)).*(sin(x(2))/x(2));
end

