clc;clear;close all
load ColorMap_me.mat

if ~exist('D:\研究生作业\数模作业2-3\线性HN','dir')
    mkdir('D:\研究生作业\数模作业2-3\线性HN')
end

% 动图参数
Fname='D:\研究生作业\数模作业2-3\线性HN\';
GifName='线性HN.gif';
delay=0.3;  % 延迟时间

f=1.e-3;  % 科氏参数
bta=1.e-7;  % 埃克曼摩擦
g=9.8;

L=1.e6;  % 边长
Nx=100;
Ny=Nx;
dx=L/Nx;
xc=0.5*Nx;

% 初值
u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
u1=u;
v1=v;

h=zeros(Nx,Ny);
[X,Y]=meshgrid(1:Nx,1:Ny);
H0=600;  % 假设平均水深
h=H0-500*exp( ( -( X-xc*1.6 ).^2-( Y-xc*1.6 ).^2 )/15^2 );
figure;
imagesc(h);
colorbar;
title('bathmetry/ m')
print([Fname,'0'],'-dpng');

V=imread([Fname,num2str(0),'.png']);
[X1,map]=rgb2ind(V,256);
imwrite(X1,map,[Fname,GifName],'gif','LoopCount',inf,'DelayTime',delay)

dt=bta/f^2*0.5;  % 时间步长
Coe_uv=1./(1+bta*dt);  % 系数

eta=zeros(Nx,Ny);
eta1=eta;
eta=0.2*exp( (-(X-xc).^2-(Y-xc).^2)/4^2 );  % 初始海表高度

ts=sum(eta(:));  % 用于检查质量守恒
disp(['initial mass is ',num2str(ts),'!']);

niter=0;
while niter<400
    j=2:Nx-1;  % 向东
    k=2:Ny;  % 向南
    u1(k,j)=u(k,j)+f*0.25*dt*(v(k-1,j)+v(k-1,j-1)+v(k,j)+v(k,j-1))-g*dt/dx*(eta(k,j)-eta(k,j-1));
    u1(k,j)=u1(k,j)*Coe_uv;
    
    j=1:Nx-1;
    k=2:Ny-1;
    v1(k,j)=v(k,j)-f*0.25*dt*(u(k,j)+u(k,j+1)+u(k+1,j)+u(k+1,j+1))-g*dt/dx*(eta(k,j)-eta(k+1,j));
    v1(k,j)=v1(k,j)*Coe_uv;
     
    j=1:Nx-1;
    k=2:Ny;
    
    eta1(k,j)=eta(k,j)-dt/dx*(0.5*(h(k-1,j+1)+h(k,j+1)).*u1(k,j+1)- 0.5*(h(k-1,j)+h(k,j)).*u1(k,j)+...
        0.5*(h(k-1,j)+h(k-1,j+1)).*v1(k-1,j)-  0.5*(h(k,j)+h(k,j+1)).*v1(k,j));
    u=u1;
    v=v1;
    
    j=2:Nx-1;
    k=2:Ny;
    u1(k,j)=u1(k,j).*(eta(k,j)+eta(k,j-1))*0.5;
    
    j=1:Nx-1;
    k=2:Ny-1;
    v1(k,j)=v1(k,j).*(eta(k,j)+eta(k+1,j))*0.5;
    
    j=1:Nx-1;
    k=2:Ny;
    eta1(k,j)=eta1(k,j)-dt/dx*(u1(k,j+1)-u1(k,j)+v1(k-1,j)-v1(k,j));
    
    ts=[ts sum(eta1(:))];
    
    eta=eta1;
    
    drawnow
    if mod(niter,2)==0 
        imagesc(eta);
        % contourf(eta)
        title(['\eta at ',num2str(niter),' steps'])
        colorbar;
        colormap(ColorMap_me)
        % set(gca,'zlim',[-1  1]*1.e0)
        print([Fname,num2str(niter)],'-dpng');

        V=imread([Fname,num2str(niter),'.png']);
        [X1,map]=rgb2ind(V,256);
        imwrite(X1,map,[Fname,GifName],'gif','WriteMode','append','DelayTime',delay)
        pause(0.15)
        % hold on 
        % quiver(u,v)
    end
    
    niter=niter+1;

end


figure;
plot(0:niter,ts,'bo-')
title('time series of mass')
xlabel('steps')
print([Fname,'time series of mass'],'-dpng');

figure('position',[10,10,800,400])
subplot(1,2,1)
imagesc(u)
title('u')
subplot(1,2,2)
imagesc(v)
title('v')
print([Fname,'u & v'],'-dpng');