clc;clear;close all

FileHead='D:\研究生作业\数模作业4\';
if ~exist('D:\研究生作业\数模作业4','dir')
    mkdir('D:\研究生作业\数模作业4')
end

%% 积分值
n=2;
S1(1)=log(6)-log(5);
for i=1:8
    S1(n)=1/n-5*S1(n-1);
    n=n+1;
end

%% 算法
n1=1:8;
S2(n1)=(1/2).*( 1./(5*(n1+1)) + 1./(6*(n1+1)) );


%% 绘图
figure
subplot(1,2,1)
plot(S1,'DisplayName','直接积分值')
hold on
plot(S2,'DisplayName','利用上下限平均算法')
xlim([1 8])
xlabel('n')
legend
title('两种算法原值')
%%%消除白色边框%%%
ax = gca;
outerpos = ax.OuterPosition; % [0, 0, 1, 1]

ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width*0.99 ax_height];
%%%消除白色边框%%%
set(gca,'fontsize',15,'fontweight','bold')


subplot(1,2,2)
plot(S1(1:8)-S2)
xlim([1 8])
xlabel('n')
title('两种算法差值')
%%%消除白色边框%%%
ax = gca;
outerpos = ax.OuterPosition; % [0, 0, 1, 1]

ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width*0.99 ax_height];
%%%消除白色边框%%%
set(gca,'fontsize',15,'fontweight','bold')

set(gcf,'position',[0 0 1920 1080])
print(gcf,'-djpeg','-r600',[FileHead,'两种方法差值.jpg'])
