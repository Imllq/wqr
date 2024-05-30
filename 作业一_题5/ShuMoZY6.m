%% 蛙跃格式，计算对流方程，u>0
clc;clear;close all

if ~exist('D:\研究生作业\数模作业6\蛙跃','dir')
    mkdir('D:\研究生作业\数模作业6\蛙跃')
end

Fname='D:\研究生作业\数模作业6\蛙跃\';
GifName='CTCS.gif';
delay=0.3;  % 延迟时间

L=1;  % 棒子的长度
Nx=181;  % 水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1.e-2;
A=5;
dt=0.2*dx/u; 
lambda=u*dt/dx;
T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );
T1=zeros(size(T0));
T2=T1;
xi=2:Nx-1;

% leap-frog第一步启动时候，需要从0到1
T1(xi)=T0(xi)-0.5*lambda*( T0(xi+1)-T0(xi-1) );
% 边界条件
T1(1)=T0(1)-0.5*lambda*( T0(2)-T0(1) );
T1(end)=T0(end)-0.5*lambda*( T0(end)-T0(end-1) );

% 绘制并写入第一张图
figure
H=area(x,T0,-A);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[-A A])
title('0s')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

count=0;
while count<6000
    count=count+1;
    T2(xi)=T0(xi)-lambda*( T1(xi+1)-T1(xi-1) );
  
    % 边界条件
    T2(1)=T0(1)-lambda*( T1(2)-T1(end) );
    T2(end)=T0(end)-lambda*( T1(1)-T1(end-1) );
    
    T0=T1;
    T1=T2;

    if mod(count,30)==0
        area(x,T2,-A);
        % plot(x,T1)
        xlim([0 L]);
        ylim([-A A]);
        title([num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')

        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
    end
end
