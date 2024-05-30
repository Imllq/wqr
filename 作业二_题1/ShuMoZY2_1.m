clc;clear;close all

if ~exist('D:\研究生作业\数模作业2-1','dir')
    mkdir('D:\研究生作业\数模作业2-1')
end

%% 2-5-10(1)

if ~exist('D:\研究生作业\数模作业2-1\2-5-10-1\','dir')
    mkdir('D:\研究生作业\数模作业2-1\2-5-10-1\')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-1\2-5-10-1\';
GifName='non_line.gif';
delay=0.3;  % 延迟时间

dx=0.1;
dt=1e-4;
x=0:0.1:15;
delta=10;
u0=sin(2*pi*x);

figure
plot(x,u0)
% hold on
u1=zeros(size(u0));
u2=zeros(size(u0));
u1(2:end-1)=u0(2:end-1)+dt*0.5*( (u0(2:end-1) + u0(3:end) ).^2 - ( u0(2:end-1) + u0(1:end-2)).^2 )/(4*dx);

print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

drawnow;
cc=0;

%leap frog
while cc<10000
    cc=cc+1;
    u2(2:end-1)=u0(2:end-1)-dt*((u1(2:end-1)+u1(3:end)).^2-(u1(2:end-1)+u0(1:end-2)).^2)/(4*dx);
    u0=u1;
    u1=u2;
    
    if mod(cc,100)==0
        plot(x,u1);
        set(gca,'ylim',[-0.2,1],'xlim',[0 12.8])
        title([num2str(cc*dt),'s'])

        print([Fname,num2str(cc)],'-dpng');

        V=imread([Fname,num2str(cc),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
    end

end
    

%% 2-5-11(1)
close all

if ~exist('D:\研究生作业\数模作业2-1\2-5-11-1\','dir')
    mkdir('D:\研究生作业\数模作业2-1\2-5-11-1\')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-1\2-5-11-1\';
GifName='non_line.gif';
delay=0.3;  % 延迟时间

dx=0.1;
dt=0.004;
x=0:0.1:15;
delta=10;
u0=sin(2*pi*x);

figure
plot(x,u0)
hold on
u1=zeros(size(u0));
ub=zeros(size(u0));
u2=zeros(size(u0));

print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

drawnow;
cc=0;
%leap frog
while cc<1000
    cc=cc+1;
    ub=0.5*(u1+u0);
    u2(2:end-1)=u0(2:end-1)-dt*(ub(2:end-1)+ub(3:end)-ub(1:end-2))+(ub(3:end).^2-ub(1:end-2).^2)/(3*dx);
    u0=u1;
    u1=u2;
    
    if mod(cc,10)==0
        plot(x,u2);
        set(gca,'ylim',[-0.2,1],'xlim',[0 12.8])
        title([num2str(cc*dt),'s'])
        
        print([Fname,num2str(cc)],'-dpng');

        V=imread([Fname,num2str(cc),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
    end

end

%% 2-5-10(2)
close all

if ~exist('D:\研究生作业\数模作业2-1\2-5-10-2\','dir')
    mkdir('D:\研究生作业\数模作业2-1\2-5-10-2\')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-1\2-5-10-2\';
GifName='non_line.gif';
delay=0.3;  % 延迟时间

dx=0.1;
dt=0.004;
x=0:0.1:15;
delta=10;
u0=1.5+sin(2*pi*x);

figure
plot(x,u0)
hold on
u1=zeros(size(u0));
u2=zeros(size(u0));
u1(2:end-1)=u0(2:end-1)+dt*0.5*( (u0(2:end-1)+u0(3:end)).^2 - (u0(2:end-1)+u0(1:end-2)).^2  )/(4*dx);

print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)


drawnow;
cc=0;
%leap frog
while cc<1000
    cc=cc+1;
    u2(2:end-1)=u0(2:end-1)-dt*((u1(2:end-1)+u1(3:end)).^2-(u1(2:end-1)+u0(1:end-2)).^2)/(4*dx);
    u0=u1;
    u1=u2;
    
    if mod(cc,10)==0
        plot(x,u1);
        set(gca,'ylim',[0.5,2.5],'xlim',[0 12.8])
        title([num2str(cc*dt),'s'])
        
        print([Fname,num2str(cc)],'-dpng');

        V=imread([Fname,num2str(cc),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
    end

end
    
%% 2-5-11(2)
close all

if ~exist('D:\研究生作业\数模作业2-1\2-5-11-2\','dir')
    mkdir('D:\研究生作业\数模作业2-1\2-5-11-2\')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-1\2-5-11-2\';
GifName='non_line.gif';
delay=0.3;  % 延迟时间

dx=0.1;
dt=0.004;
x=0:0.1:15;
delta=10;
u0=1.5+sin(2*pi*x);

figure
plot(x,u0)
hold on
u1=zeros(size(u0));
ub=zeros(size(u0));
u2=zeros(size(u0));

print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)


drawnow;
cc=0;
%leap frog
while cc<1000
    cc=cc+1;
    ub=0.5*(u1+u0);
    u2(2:end-1)=u0(2:end-1)-dt*(ub(2:end-1)+ub(3:end)-ub(1:end-2))+(ub(3:end).^2-ub(1:end-2).^2)/(3*dx);
    u0=u1;
    u1=u2;
    
    if mod(cc,10)==0
        plot(x,u2);
        set(gca,'ylim',[0.5,2.5],'xlim',[0 12.8])
        title([num2str(cc*dt),'s'])
        
        print([Fname,num2str(cc)],'-dpng');

        V=imread([Fname,num2str(cc),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
    end

end
