clc;clear;close all
load ColorMap_me.mat

if ~exist('D:\研究生作业\数模作业2-2\有beta涡度','dir')
    mkdir('D:\研究生作业\数模作业2-2\有beta涡度')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-2\有beta涡度\';
GifName='beta.gif';
delay=0.3;  % 延迟时间

dx=1;  % 网格间距 
N=100;  % 方形网格，每维181个点
ix=2:N-1;  % 内部点
jy=2:N-1;

coeffi_Jacobian= -1/(12*dx^2);  % Arakawa 1966 Eq(45,46)

r_d2x=1/(2*dx);
r_dx2=1/(dx^2);
Nts=8000;  % 时间步数

Ah=1.e-5;  % 粘性系数
bta=1e-2;  % beta效应强度
Amp=1;

% 涡度场，3个时间层
zeta0= zeros(N,N);  % 初始场，第n-1层
zeta1= zeros(N,N);  % 第n层
zeta2= zeros(N,N);  % 第n+1层
psi=zeros(N);  % 储存泊松方程的解
psi_guess=rand(N);  % 初猜值

% 初始化原始相对涡度场
for i=1:N
    for j=1:N
        r2=(j/N-0.35)^2*60+(i/N-0.25)^2*100;
        r3=(j/N-0.65)^2*30+(i/N-0.75)^2*100;
        zeta0(i,j)=exp(-1*r2)-exp(-1*r3);
    end
end
zeta0=Amp*zeta0 ;
[xx,yy]=meshgrid(1:N,1:N);

figure;
imagesc(zeta0);
title('initial vorticity'); 
colorbar


% SOR求解初始场的psi
SOR2d
figure
imagesc(psi);
title('SOR solution of \psi for the initial \zeta field')
colorbar

u=max(max(abs(diff(psi)/dx)));
dt=(dx/u)*0.2;
             
jac
zeta1=zeta0+J0*dt;  % Euler 时间向前预测下一步
zeta0=zeta1;  
SOR2d;  
jac;
J1=J0;

tic;
count=0;
ts=0;

figure('position',[0,0,1600,600])
drawnow

T=0;

subplot(1,2,1);
contourf(zeta1,50,'LineStyle','none');
colormap(ColorMap_me)
colorbar;
clim([-1 1]*Amp*0.8);
set(gca,'ydir','reverse')
title([num2str(ts),' steps, cpu time ',...
    num2str(toc,'%3.2f'),'s, model T=',num2str(T,'%3.2f'),'s']);

subplot(1,2,2)
plot(0,'bx-')
xlabel('time step')
ylabel('enstrophy')
pause(0.1);
u=max(max(abs(diff(psi)/dx)));
dt=min(dx/u)*0.1;

print([Fname,'0'],'-dpng');
V=imread([Fname,num2str(0),'.png']);
[X,map]=rgb2ind(V,256);
imwrite(X,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

while(ts<Nts)
    ts=ts+1;
    T=T+dt;

    zeta2=zeta1+dt*(1.5*J1-0.5*J0) ;  % 非线性项的Adam-Bath二阶时间格式
    
    % 有beta
    zeta2(jy,ix)=zeta2(jy,ix)+(zeta2(jy+1,ix)+zeta2(jy-1,ix)...
        +zeta2(jy,ix+1)+zeta2(jy,ix-1))*Ah*r_dx2*dt-...
        dt*bta*(psi(jy,ix+1)-psi(jy,ix-1))*0.5;

    % 无beta
    % zeta2(jy,ix)=zeta2(jy,ix)+...
    %     (zeta2(jy+1,ix)+zeta2(jy-1,ix)+zeta2(jy,ix+1)+zeta2(jy,ix-1))*Ah*r_dx2*dt; 

    diffusion_periobc

    zeta2=zeta2/(1+4*Ah*dt*r_dx2);  
    
    zeta1=zeta2;
    zeta0=zeta2;

    SOR2d;  
    jac;        
    J2=J0;
    
    J0=J1;
    J1=J2;
    
    diag_ts(ts,1)=sum(sum(zeta2.^2));
    
    cputime_per_step=toc;
    if mod(ts,20)==0
        subplot(1,2,1);
        contourf(zeta2,50,'LineStyle','none');
        colormap(ColorMap_me)
        colorbar;
        clim([-1 1]*Amp*0.8);
        set(gca,'ydir','reverse')
        title([num2str(ts),' steps, cpu time ',...
            num2str(cputime_per_step,'%3.2f'),'s, model T=',num2str(T,'%3.2f'),'s']);
   
        subplot(1,2,2)
        plot(diag_ts,'bx-')
        xlabel('time step')
        ylabel('enstrophy')
        pause(0.1);
        u=max(max(abs(diff(psi)/dx)));
        dt=min(dx/u)*0.1;

        print([Fname,num2str(ts)],'-dpng');
        V=imread([Fname,num2str(ts),'.png']);
        [X,map]=rgb2ind(V,256);
        imwrite(X,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        
    end
end


