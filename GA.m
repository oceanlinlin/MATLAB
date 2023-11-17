%------标准遗传算法求函数极值、最大值、最小值------
%------给定初始化条件-----------
clear all;
close all;
clc;
NP = 50; %群体粒子个数
L = 20;  %二进制编码长度
G = 100; %最大遗传代数
Pc = 0.8; %交叉概率
Pm = 0.1; %变异概率
Xs = 1; %上限
Xx = -1;  %下限
f = randi([0,1],NP,L); %随机初始化种群，由于染色体为二进制编码。

format long;

%--------循环迭代次数------------
for k=1:G
    %---------将二进制编码解码为定义域范围内的十进制-----------
    for i=1:NP
       U = f(i,:) ;
       m = 0 ;
       for j=1:L
         m = U(j)*2^(j-1)+m; %将二进制转化为十进制
       end
       x(i) = Xx + m*(Xs-Xx)/(2^L-1); %矩阵x用于存储十进制的数值，解码到定义域中
       Fit(i) = func3(x(i));  %初始时每个样本的适应度值，等于目标函数的值
    end
    maxFit = max(Fit);  %最大值
    minFit = min(Fit);  %最小值
    rr = find(Fit == maxFit);  %返回最大/小适应度值所在位置,返回的是一个数组
    fbest = f(rr(1,1),:); %返回最优个体的初始种群对应的值
    xbest = x(rr(1,1)); %返回最优适应度对应的染色体十进制的值
    Fit = (Fit - minFit)/(maxFit - minFit); %归一化适应度值

    %---------基于轮盘赌的复制操作-----------
    sum_Fit = sum(Fit); %计算种群中所有种群的适应度值
    fitvalue = Fit./sum_Fit; %每个种群被选中的概率
    fitvalue = cumsum(fitvalue); %计算每个种群的累计选择概率
    ms = sort(rand(NP,1)); %随机生成(0,1)内的数值，将该随机数作为选择指针来确定被选个体，并将其升序排列
    fiti = 1;
    newi = 1;
    while newi<=NP
        if(ms(newi)<fitvalue(fiti)) %随机新种群的概率小于种群选择概率
            nf(newi,:) = f(fiti,:); %新种群的第newi行 = 此次选择的种群
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end

    %---------基于概率的交叉操作-----------
    for i = 1:NP-1
        p = rand; %随机生成一个[0,1]的概率p
        if p<Pc %控制交叉的颜色体数
            q = randi([0,1],1,L); %随机生成要交叉的基因位置
            for j = 1:L
                if q(j) == 1 %对应的值为1，找到该点所在的位置视为交叉位置
                    temp = nf(i+1,j);
                    nf(i+1,j) = nf(i,j); %基因互换
                    nf(i,j) = temp;
                end
            end
        end
    end

    %---------基于概率的变异操作-----------
    i = 1;
    while i<= round(NP*Pm) %round为四舍五入的整数，NP*Pm表示NP个种群中以Pm为概率变异的种群数
        h = randi([1,NP],1,1); %生成一个1*1的矩阵，矩阵元素在[1,NP]之间，即随机选取一个需要变异的染色体
        for j = 1:round(L*Pm)
            g = randi([1,L],1,1); %随机选取需要变异的基因数
            nf(h,g) =~nf(h,g); %变异基因取反
        end
        i = i+1;
    end
    f = nf;  %输出新一代种群
    f(1,:) = fbest; %保留最优个体在新种群中
    trace1(k) = maxFit; %历代最优适应度值
    %trace2(k) = mean(Fit); %平均适应度值
    %trace3(k) = minFit; %最小值
end

xbest; %最优个体
figure;
%hold on
plot(1:k,trace1(1:k),'r')
%plot(1:k,trace2(1:k),'b')
%plot(1:k,trace3(1:k),'g')
title('最优个体适应度','fontsize',12);
xlabel('进化代数','FontSize',12);
ylabel('适应度','FontSize',12);
%legend('最大值','平均适应值','最小值');

%hold off

%以下是保存为GIF的程序
%frame = getframe(gcf);
%imind = frame2im(frame);
%[imind,cm] = rgb2ind(imind,256);
%if i == 1
%    imwrite(imind,cm,'Z1.gif','gif','LoopCount',inf,'DelayTime',0);
%else
%    imwrite(imind,cm,'Z1.gif','gif','WriteMode','overwrite','DelayTime',0);
%end


%------最优个体-----
%fprintf('最优个体编码：\n')
%disp(fbest)

%------每次迭代的最优值-----
%fprintf('每次迭代的最优值：\n')
%disp(trace1)

%------每次迭代的平均适应度值----
%fprintf('每次迭代的平均适应度值:\n')
%disp(trace2)

%------每次迭代的最差值----
%fprintf('每次迭代的最差值:\n')
%disp(trace3)

 %-----------适应度函数-------------
function result = func3(x)
result = 10+x.*cos(5*pi*x);
end



