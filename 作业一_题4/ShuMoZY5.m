%% BTCS-dt
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.E-6;  % 热传导系数
A=1;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业5\BTCS','dir')
    mkdir('D:\研究生作业\数模作业5\BTCS')
end

% 动图参数
Fname='D:\研究生作业\数模作业5\BTCS\';
GifName='BTCS-dt.gif';
delay=0.3;  % 延迟时间

% 初始时刻，t=0时温度分布
% T0=A*cos((pi/L).*x);  % 按照初始条件得出T0温度分布
% T0=0.5*A*(1+cos(1*x/L*pi));  % PPT中T0温度分布
% T0=A*exp( - (x-0.5*L).^4*1.e5 );
% T0=A*cos(pi/L.*x);
T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );

D=zeros(Nx-2,Nx-2);
for i=1:Nx-2
    D(i,i)=-(1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu;
    end
    if i>1
        D(i,i-1)=mu;
    end
end

% 第一类边界条件：T(0)=0=T(N)

% 第二类边界条件：在隐式数组形式上T_x=0
D(1,1)=-(1+mu);
D(end,end)=-(1+mu);

D_inv= D^-1;

% 绘制并写入第一张图
figure
H=area(x,T0,0);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第0步-BTCS-dt')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值
% T2=zeros(size(T0));  % 储存t=2时刻温度计算值


% 开始计算
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while count<20000
    count=count+1;

    % T1(xi)=(1-2*mu)*T0(xi) + mu*( T0(xi+1)+T0(xi-1) );
    T1(xi)= - D_inv*T0(xi)';
    T1(1)=T0(1)+dt/dx^2*K*T1(2) ;

    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;

    if mod(count,10)==0

        area(x,T1,0)
        set(gca,'xlim',[0 L],'ylim',[0 A])
        title(['第',num2str(count),'步-BTCS-dt'])
        xlabel('x')
        ylabel('T')
        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
        
    end
    T1=zeros(size(T0));

end


%% FTCS-dt
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.E-6;  % 热传导系数
A=1;  % 振幅
dt=10;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业5\FTCS','dir')
    mkdir('D:\研究生作业\数模作业5\FTCS')
end

% 动图参数
Fname='D:\研究生作业\数模作业5\FTCS\';
GifName='FTCS-dt.gif';
delay=0.3;  % 延迟时间

% 初始时刻，t=0时温度分布
% T0=A*cos((pi/L).*x);  % 按照初始条件得出T0温度分布
% T0=0.5*A*(1+cos(1*x/L*pi));  % PPT中T0温度分布
% T0=A*exp( - (x-0.5*L).^4*1.e5 );
% T0=A*cos(pi/L.*x);
T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );

% 绘制并写入第一张图
figure
H=area(x,T0,0);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第0步-FTCS-dt')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值
% T2=zeros(size(T0));  % 储存t=2时刻温度计算值


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

        area(x,T1,0)
        % refreshdata(H, 'caller');
        % hold on
        % drawnow
        set(gca,'xlim',[0 L],'ylim',[0 A])
        title(['第',num2str(count),'步-BTCS-dt'])
        xlabel('x')
        ylabel('T')
        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
        
    end
    T1=zeros(size(T0));

end


%% BTCS-2dt
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.E-6;  % 热传导系数
A=1;  % 振幅
dt=20;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业5\BTCS','dir')
    mkdir('D:\研究生作业\数模作业5\BTCS')
end

% 动图参数
Fname='D:\研究生作业\数模作业5\BTCS\';
GifName='BTCS-2dt.gif';
delay=0.3;  % 延迟时间

% 初始时刻，t=0时温度分布
% T0=A*cos((pi/L).*x);  % 按照初始条件得出T0温度分布
% T0=0.5*A*(1+cos(1*x/L*pi));  % PPT中T0温度分布
% T0=A*exp( - (x-0.5*L).^4*1.e5 );
% T0=A*cos(pi/L.*x);
T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );

D=zeros(Nx-2,Nx-2);
for i=1:Nx-2
    D(i,i)=-(1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu;
    end
    if i>1
        D(i,i-1)=mu;
    end
end

% 第一类边界条件：T(0)=0=T(N)

% 第二类边界条件：在隐式数组形式上T_x=0
D(1,1)=-(1+mu);
D(end,end)=-(1+mu);

D_inv= D^-1;

% 绘制并写入第一张图
figure
H=area(x,T0,0);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第0步-BTCS-2dt')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值
% T2=zeros(size(T0));  % 储存t=2时刻温度计算值


% 开始计算
xi=2:Nx-1;  % 头尾格点由初始条件和边界条件给出，xi为需计算格点
while count<20000
    count=count+1;

    % T1(xi)=(1-2*mu)*T0(xi) + mu*( T0(xi+1)+T0(xi-1) );
    T1(xi)= - D_inv*T0(xi)';
    T1(1)=T0(1)+dt/dx^2*K*T1(2) ;

    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;

    if mod(count,10)==0

        area(x,T1,0)
        set(gca,'xlim',[0 L],'ylim',[0 A])
        title(['第',num2str(count),'步-BTCS-2dt'])
        xlabel('x')
        ylabel('T')
        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
        
    end
    T1=zeros(size(T0));

end


%% FTCS-2dt
clc;clear;close all

% 输入初始条件
% 预设条件
L=1;  % 棒子长度
Nx=181;  % 水平节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);  % 格点间距
K=1.E-6;  % 热传导系数
A=1;  % 振幅
dt=20;  % 时间步长
mu=(K*dt)/(dx^2);

if ~exist('D:\研究生作业\数模作业5\FTCS','dir')
    mkdir('D:\研究生作业\数模作业5\FTCS')
end

% 动图参数
Fname='D:\研究生作业\数模作业5\FTCS\';
GifName='FTCS-2dt.gif';
delay=0.3;  % 延迟时间

% 初始时刻，t=0时温度分布
% T0=A*cos((pi/L).*x);  % 按照初始条件得出T0温度分布
% T0=0.5*A*(1+cos(1*x/L*pi));  % PPT中T0温度分布
% T0=A*exp( - (x-0.5*L).^4*1.e5 );
% T0=A*cos(pi/L.*x);
T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );

% 绘制并写入第一张图
figure
H=area(x,T0,0);
grid minor
drawnow
set(gca,'xlim',[0 L],'ylim',[0 A])
title('第0步-FTCS-2dt')
xlabel('x')
ylabel('T')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

% 初始参数
count=0;

T1=zeros(size(T0));  % 储存t=1时刻温度计算值
% T2=zeros(size(T0));  % 储存t=2时刻温度计算值


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

        area(x,T1,0)
        % refreshdata(H, 'caller');
        % hold on
        % drawnow
        set(gca,'xlim',[0 L],'ylim',[-A A])
        title(['第',num2str(count),'步-BTCS-2dt'])
        xlabel('x')
        ylabel('T')
        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
        
    end
    T1=zeros(size(T0));

end

