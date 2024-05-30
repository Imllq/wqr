clc;clear;close all
load ColorMap_me.mat

if ~exist('D:\研究生作业\数模作业2-4\线性HN开','dir')
    mkdir('D:\研究生作业\数模作业2-4\线性HN开')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-4\线性HN开\';
GifName='线性HN开.gif';
delay=0.3;  % 延迟时间

f=1.e-4;  % 科氏参数
bta=1.e-7;  % 埃克曼摩擦
g=9.8;
L=1.e6;  % 边长
Nx=200;  % 分辨率
Ny=Nx;
dx=L/Nx;  % 间距
xc=0.5*Nx;  % 初始扰动位置

u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
u1=u;           
v1=v;          
eta=zeros(Nx,Ny); 
eta1=eta;         
h=zeros(Nx,Ny);   

[X,Y]=meshgrid(1:Nx,1:Ny);
H0=2000;
h=ones(Nx,Ny)*H0;

eta=0.05*exp( -(X-xc).^2/4^2-(Y-xc).^2/8^2);

dt=min([dx/sqrt(g*max(h(:))),bta/f^2]);
ns=max([10/dt, 10]);
Coe_uv=1./(1+bta*dt);

figure('Position',[0 0 1960 1400])
mesh((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,eta);
drawnow
set(gca,'zlim',[-0.02  0.05]);
colorbar
colormap(ColorMap_me)
zlabel('\eta, m')
title([num2str(round(0*dt)),' s'])
xlabel('x, km')
ylabel('y, km')
zlabel('\eta, m')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X1,map]=rgb2ind(V,256);
imwrite(X1,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

mu=sqrt(9.8*H0)*dt/dx;

ts=sum(eta(:)) ;
disp(['initial mass is ',num2str(ts),'!']);

niter=0;T=0;
while niter*dt<=3*3600
    
    j=2:Nx-1; % 向东
    k=1:Ny-1; %向北
    u1(k,j)=u(k,j)+f*0.25*dt*(v(k,j)+v(k+1,j)+v(k,j-1)+v(k+1,j-1))-g*dt/dx*(eta(k,j)-eta(k,j-1));
    u1(k,j)=u1(k,j)*Coe_uv;
    
    j=Nx;    
    u1(k,j)=(1-mu).*u(k,j)+mu.*u(k,j-1);

    j=1:Nx-1;
    k=2:Ny-1;
    v1(k,j)=v(k,j)-f*0.25*dt*(u(k,j)+u(k,j+1)+u(k-1,j)+u(k-1,j+1))-g*dt/dx*(eta(k,j)-eta(k-1,j));
    v1(k,j)=v1(k,j)*Coe_uv;
    
    j=Nx;
    v1(k,j)=(1-mu).*v(k,j)+mu.*v(k,j-1);
    
    j=1:Nx-1;
    k=1:Ny-1;
    eta1(k,j)=eta(k,j)-dt/dx*( ......
        0.5*(h(k+1,j+1)+h(k,j+1)).*u1(k,j+1)-0.5*(h(k+1,j)+h(k,j)).*u1(k,j) +.....
        0.5*(h(k+1,j)+h(k+1,j+1)).*v1(k+1,j)-0.5*(h(k,j)+h(k,j+1)).*v1(k,j));
    
    j=Nx;
    eta1(k,j)=(1-mu).*eta(k,j)+mu.*eta(k,j-1);
    
    u=u1;
    v=v1;
    eta=eta1;
   
    ts=[ts sum(eta1(:))];
    
    drawnow
    if mod(niter,ns)==0
        mesh((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,eta);
        drawnow
        set(gca,'zlim',[-0.02  0.05]);
        colorbar
        colormap(ColorMap_me)
        zlabel('\eta, m')
        title([num2str(round(niter*dt)),' s'])
        xlabel('x, km')
        ylabel('y, km')
        zlabel('\eta, m')
        print([Fname,num2str(niter)],'-dpng');

        V=imread([Fname,num2str(niter),'.png']);
        [X1,map]=rgb2ind(V,256);
        imwrite(X1,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.3)
        
    end
    
    T=T+dt;
    niter=niter+1;
end

figure;
plot( ts,'b-')
title('time series of mass')
xlabel('steps')

figure('position',[10,10,800,400])
subplot(121)
imagesc(u)
title('u')
subplot(122)
imagesc(v)
title('v')
