function [IPSO_pbestloc,IPSO_gbestloc,IPSO_Best,IPSO_end] = IPSO(N,D,T,Xmax,Xmin,Vmax,Vmin,gb)

%调整惯性权重初始化的参数
a = 1;
b = 0.7;
d = 2;

%调整学习因子初始化的参数
e1 = 1;
e2 = 0.7;
f1 = 1;
f2 = 0.7;

%------初始化种群的个体（可以在此处限定位置和速度的范围）------------
x = rand(N,D) * (Xmax-Xmin) + Xmin;
v = rand(N,D) * (Vmax-Vmin) + Vmin;


%------初始化个体最优位置和最优值----------------------

IPSO_pbestloc = x;  %因为是第一代，所以初始时个体最优位置p即为初始位置
IPSO_pbest = ones(N,1);  %构建N*1矩阵，用于存储个体最优值
for i=1:N
    IPSO_pbest(i) = func2(x(i,:));  %各个粒子的最优适应度（最优值），因为是第一代，个体最优值即为函数本身的最优适应度
end

%------初始化全局最优位置和全局最优值------------
IPSO_gbestloc = ones(1,D); %构建1*N矩阵，由于存储全局最优位置
IPSO_Best = inf;  %初始时，全局最优值设为0，之后比较时选取大于号用于求取函数最大值
%若初始时候gbest取值为inf,之后比较时用小于号来求取最小值
%论文中f1函数求最大值时候可以设置PSO_Best为0;f2函数得设置其值初始为负无穷，或者更直接一点，求最大值时候全设置为-inf
%论文中f3函数有最小值，可设置初始全局最优值为inf，下面对应的＞要改为＜
for i = 1:N
    if(IPSO_pbest(i)<IPSO_Best)  %大于符号用于求取最大值，小于符号用于求取最小值
        IPSO_Best = IPSO_pbest(i);  %此时的全局最优值为第i个粒子的最优值
        IPSO_gbestloc = IPSO_pbestloc(i,:);    %此时的全局最优位置为粒子矩阵p中第i个粒子对应的初始位置
    end
end
gb = ones(1,T);   %用于后续绘图，显示每次迭代的全局最优值


%----------------按照公式迭代直到满足精度或迭代次数----------------------
for i=1:T
        %动态调整惯性权重
        w = a+(a-b)*((i)^d)/((T)^d);
        %动态调整学习因子
        c1 = (e1-e2)*i/T + e2;
        c2 = (f1-f2)*i/T + f2;
    for j=1:N
        %-----------更新个体最优位置和最优值-------------
        if(func2(x(j,:))<IPSO_pbest(j)) %若此时的适应函数值大于之前的最优值，目的是求最大值；若为小于符号，求的是最小值
            IPSO_pbest(j) = func2(x(j,:));
            IPSO_pbestloc(j,:) = x(j,:);   %此时的个体最优位置等于取当前最优值时候的粒子位置
        end
        %-----------更新全局最优位置和最优值-------------
        if(IPSO_pbest(j)<IPSO_Best) %若为大于符号求最大值，小于符号求最小值
            IPSO_Best = IPSO_pbest(j); %此时的全局最优值等于当前粒子的个体最优值
            IPSO_gbestloc = IPSO_pbestloc(j,:);
        end
        %-----------更新位置和速度值-------------
        v(j,:) = w * v(j,:) + c1 * rand *(IPSO_pbestloc(j,:)-x(j,:)) + c2 * rand *(IPSO_gbestloc - x(j,:));
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
    gb(i) = IPSO_Best;
end
%IPSO_gbestloc = vectorToDecimal(IPSO_gbestloc)  %显示全局最优的粒子位置
IPSO_gbestloc;
IPSO_Best = gb;  %显示全局最优的粒子对应的最优值
IPSO_end = gb(end);
IPSO_pbestloc;

%hf1 = figure;
%figure;
%plot(gb,"r-.",'LineWidth',2)
%xlabel('迭代次数');
%ylabel("适应度值");
%title("适应度进化曲线");
%legend('IPSO');
%savefig(hf1,'IPSO曲线.fig');


end

function result  = func2(x)
result = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
%result = x(1)^2 + x(2)^2;
%result = 20+(x(1)^2-10*cos(2*pi*x(1)))+(x(2)^2-10*cos(2*pi*x(2)));
%result = 20*exp(-0.2*sqrt((x(1)^2+x(2)^2)/2)) + exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2) - exp(1) - 20;
%result = (sin(x(1))/x(1)).*(sin(x(2))/x(2));
end


