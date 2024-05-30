%% 数值解
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.E-6;  % 热传导系数
lambda=pi/L;
A=1;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业3\FTCS','dir')
    mkdir('D:\研究生作业\数模作业3\FTCS')
end

% 动图参数
Fname='D:\研究生作业\数模作业3\FTCS\';
GifName='FTCS.gif';
delay=0.3;  % 延迟时间

% 设置T0，并计算
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );
an=2/Nx*fft(T0)';
kn=[0:Nx/2  , -Nx/2+1:-1]'.*2*pi/L; %L=Nx*dx
expikx= exp(1j.*(kn*x));

% 绘制并写入第一张图
figure
H=area(x,T0,-A);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[-A A])
title('第0步')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值

% 开始计算
% 使用FTCS(Forward Time Central Space)
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while count<20000
    count=count+1;

    % T1(xi)=(1-2*mu)*T0(xi) + mu*( T0(xi+1)+T0(xi-1) );
    T1(xi)=T0(xi) + mu.*( T0(xi+1)-2.*T0(xi)+T0(xi-1) );

    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;

    if mod(count,200)==0
        % figure(1)
        % grid minor
        % hold on

        area(x,T1,-A)
        % refreshdata(H, 'caller');
        % hold on
        % drawnow
        set(gca,'xlim',[0 L],'ylim',[-A A])
        title(['第',num2str(count),'步'])
        xlabel('x')
        ylabel('T')
        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
    T1=zeros(size(T0));

end


%% 解析解
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,1,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.e-6;  % 热传导系数
lambda=pi/L;
A=1;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业3\JX','dir')
    mkdir('D:\研究生作业\数模作业3\JX')
end

% 动图参数
Fname='D:\研究生作业\数模作业3\JX\';
GifName='JieXi.gif';
delay=0.3;  % 延迟时间

% 初始时刻，t=0时温度分布
t=0;
T0=A*exp(-(x-0.5*L).^2/(0.1*L)^2);
an=2/Nx*fft(T0)';
kn=[0:Nx/2,-Nx/2+1:-1]'.*2*pi/L;
kn(181)=0;
expikx=exp(1j.*(kn*x));

% 绘制并写入第一张图
figure
T_exact=( an.*exp(-K*kn.^2*0*dt) )'*expikx;
area(x,real(T_exact),-A)
grid minor
% drawnow
set(gca,'xlim',[0 L],'ylim',[-A A])
title('第0步')
xlabel('x')
ylabel('T')
print([Fname,'J-0'],'-dpng');

V=imread([Fname,'J-',num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 开始计算
% 使用数值解
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while t<200000
    t=t+dt;

    T_exact=( an.*exp(-K*kn.^2*t*dt) )'*expikx;

    if mod(t,2000)==0
        % figure(1)
        grid minor
        % hold on

        area(x,real(T_exact),-A)
        % hold on
        drawnow
        set(gca,'xlim',[0 L],'ylim',[-A A])
        title(['第',num2str(t/10),'步'])
        xlabel('x')
        ylabel('T')
        print([Fname,'J-',num2str(t/10)],'-dpng');

        V=imread([Fname,'J-',num2str(t/10),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)

        pause(0.15)
        
    end

end

