%% 迎风格式，计算对流方程，u>0
clc;clear;close all

if ~exist('D:\研究生作业\数模作业8\迎风','dir')
    mkdir('D:\研究生作业\数模作业8\迎风')
end

Fname='D:\研究生作业\数模作业8\迎风\';
GifName='FTCS.gif';
delay=0.3;  % 延迟时间

L=1;  % 棒子的长度
Nx=181;  % 水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1.e-2;  % 1cm/s
A=5;
dt=0.1*dx/u;  % 适当扩大观察变化
lambda=u*dt/dx;
K=1.E-6;
mu=(K*dt)/(dx^2);

T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );
figure
area(x,T0,-A);
T1=zeros(size(T0));
xi=2:Nx-1;

drawnow;
set(gca,'xlim',[0 L],'ylim',[-A A])
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

count=0;
while count<3000
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+(mu-lambda/2)*T0(xi+1)+(mu+lambda/2)*T0(xi-1);
  
    % 边界条件
    T1(1)=T1(2);
    T1(end)=T1(end-1);

    T0=T1;

    if mod(count,30)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title([num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')

        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
end


%% BTCS格式，计算对流方程，u>0
clc;clear;close all

if ~exist('D:\研究生作业\数模作业8\BTCS','dir')
    mkdir('D:\研究生作业\数模作业8\BTCS')
end

Fname='D:\研究生作业\数模作业8\BTCS\';
GifName='BTCS.gif';
delay=0.3;  % 延迟时间

L=1;  % 棒子的长度
Nx=181;  % 水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1.e-2;  % 1cm/s
A=5;
dt=0.1*dx/u;  % 适当扩大观察变化
lambda=u*dt/dx;
K=1.E-6;
mu=(K*dt)/(dx^2);

T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );
figure
area(x,T0,-A);
T1=zeros(size(T0));
xi=2:Nx-1;

drawnow;
set(gca,'xlim',[0 L],'ylim',[-A A])
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

count=0;
while count<3000
    count=count+1;

    D=zeros(Nx-2,Nx-2);
    for i=1:Nx-2
        D(i,i)= - (1+2*mu);
        if i<Nx-2
            D(i,i+1)=mu-lambda/2;
        end
        if i>1
            D(i,i-1)=mu+lambda/2;
        end
    end

    D(1,1)=-(1+mu);
    D(end,end)=-(1+mu);

    D_inv= D^-1;

    % T1(1)=T0(1)+dt/dx^2*K*T1(2) ;
  
    % 边界条件
    T1(1)=T1(2);
    T1(end)=T1(end-1);

    T0=T1;

    if mod(count,30)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title([num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')

        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
end


%% 迎风格式，计算对流方程，u>0
clc;clear;close all

if ~exist('D:\研究生作业\数模作业8\迎风','dir')
    mkdir('D:\研究生作业\数模作业8\迎风')
end

Fname='D:\研究生作业\数模作业8\迎风\';
GifName='FTCS.gif';
delay=0.3;  % 延迟时间

L=1;  % 棒子的长度
Nx=181;  % 水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1.e-2;  % 1cm/s
A=5;
dt=0.1*dx/u;  % 适当扩大观察变化
lambda=u*dt/dx;
% while lambda<=1
%     dt=dt+0.1;
%     lambda=u*dt/dx;
% end
K=1.E-6;
mu=(K*dt)/(dx^2);

T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );
figure
area(x,T0,-A);
T1=zeros(size(T0));
xi=2:Nx-1;

drawnow;
set(gca,'xlim',[0 L],'ylim',[-A A])
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

count=0;
while count<3000
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+(mu-lambda/2)*T0(xi+1)+(mu+lambda/2)*T0(xi-1);
  
    % 边界条件
    T1(1)=T1(2);
    T1(end)=T1(end-1);

    T0=T1;

    if mod(count,30)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title([num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')

        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
end


%% BTCS格式，计算对流方程，u>0
clc;clear;close all

if ~exist('D:\研究生作业\数模作业8\BTCS','dir')
    mkdir('D:\研究生作业\数模作业8\BTCS')
end

Fname='D:\研究生作业\数模作业8\BTCS\';
GifName='BTCS.gif';
delay=0.3;  % 延迟时间

L=1;  % 棒子的长度
Nx=181;  % 水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1.e-2;  % 1cm/s
A=5;
dt=0.1*dx/u;  % 适当扩大观察变化
lambda=u*dt/dx;
while lambda<=1
    dt=dt+0.1;
    lambda=u*dt/dx;
end
K=1.E-6;
mu=(K*dt)/(dx^2);

T0=A*exp( -(x-0.5*L).^2/(0.1*L)^2 );
figure
area(x,T0,-A);
T1=zeros(size(T0));
xi=2:Nx-1;

drawnow;
set(gca,'xlim',[0 L],'ylim',[-A A])
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

count=0;
while count<3000
    count=count+1;

    D=zeros(Nx-2,Nx-2);
    for i=1:Nx-2
        D(i,i)= - (1+2*mu);
        if i<Nx-2
            D(i,i+1)=mu-lambda/2;
        end
        if i>1
            D(i,i-1)=mu+lambda/2;
        end
    end

    D(1,1)=-(1+mu);
    D(end,end)=-(1+mu);

    D_inv= D^-1;

    % T1(1)=T0(1)+dt/dx^2*K*T1(2) ;
  
    % 边界条件
    T1(1)=T1(2);
    T1(end)=T1(end-1);

    T0=T1;

    if mod(count,30)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title([num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')

        print([Fname,num2str(count)],'-dpng');

        V=imread([Fname,num2str(count),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
end


