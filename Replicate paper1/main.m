%-----主函数-----
clear;
clc;
close all;

%-----初始化参数-----
N = 20; %群体粒子个数
D = 10;  %粒子搜索空间维数
T = 500; %最大迭代次数
Xmax = 100; %位置最大值
Xmin = -100; %位置最小值
Vmax = 10;  %速度最大值
Vmin = -10; %速度最小值


k = 20;%试验次数
for i=1:k
    [PSO_pbestloc,PSO_gbestloc,PSO_Best,PSO_end] = PSO(N,D,T,Xmax,Xmin,Vmax,Vmin);
    [IPSO_pbestloc,IPSO_gbestloc,IPSO_Best,IPSO_end] = IPSO(N,D,T,Xmax,Xmin,Vmax,Vmin);
    %PSO_average(i,:) = PSO_gbestloc;
    %IPSO_avgerage(i,:) = IPSO_gbestloc;
    A(i,:)=PSO_gbestloc;%第i次实验全局最优位置矩阵
    B(i,:)=IPSO_gbestloc;
    C(i)=PSO_end;%第i次实验全局最优值
    E(i) = IPSO_end;
 
end




%plot(PSO_Best,'-b','LineWidth',2)
semilogy(PSO_Best,'--k','LineWidth',2)
hold on
%plot(IPSO_Best,'-.r','LineWidth',2)
semilogy(IPSO_Best,'-.r','LineWidth',2)
hold on
legend('PSO','IPSO');
xlabel("迭代次数");
ylabel("适应度值");
title("适应度曲线");

