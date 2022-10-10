function plot_rain3D

nx = 128; ny=128;  m=100+1;
%dt=0.01;

Us = 100d0/9d0; 
Ts = 15d0/60d0; 
Ths = 3.0d0; 
Ls = 10d0;
qs = 10d0; 
    
movies=12*60+1; %100 plus the one at t=0

epsilon=0.1;

Tfinal=(1/60)*(movies-1);
%T=0:dt:Tfinal;

DeltaT=1/60;

Lx = 128;
Ly = 128;
Lz = 15;

dx = Lx/nx; dy=Ly/ny; dz=Lz/(m-1);

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy; %dy/2:dy:Ly-dy/2;
z = 0:dz:Lz ;

z_w = 0:dz:Lz ;
z_u = dz/2:dz:Lz+dz/2 ;


[xGrid, yGrid, zGrid] = meshgrid(y,x,z);

a=10;
ubg=zeros(m,1);
vbg=zeros(m,1);
for k=1:m
    zk=(k-0.5)*dz;
    
    if zk< 12
        ubg(k)=a*(cos(pi*zk/12)-cos(2*pi*zk/12));
    else
        ubg(k)=-2*a;
    end
    
    if zk >= 11d0
		vbg(k) = 8;
    elseif zk >= 5
		vbg(k)=8*(zk-5)/(11-5);
    end
end

% plot(vbg,z)
% stop

%Approximate velocity for squall line
u_squall=max(ubg); %1.125*a;
u_squall=u_squall*40; %In km/hr

%NoSteps=Tfinal*100+1;
NoStepsMov=Tfinal/(movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1);
end



%%

% qv = load('Data2D/QvYHalf_99.dat');
% movies = size(qv,1)/(nx*m)
% qv = reshape(qv,nx,m,movies); %This are already in m/s

qr = load('Data2D/QrYHalf_99.dat');
movies = size(qr,1)/(nx*m)
qr = reshape(qr,nx,m,movies); %This are already in m/s

qt = qr; %qv*0+qr;

MovAvi = VideoWriter('Rain.avi');
MovAvi.Quality = 100;
MovAvi.FrameRate = 5;
open(MovAvi)

figure(1)
hold on
for mov=1:1:60*1 %movies
    figure(1)
    clf
    contour(x,z,transpose(qt(:,:,mov)),100); colorbar;
    caxis([.1 5])
    print('Frame','-depsc',figure(1))
    
    F = getframe(figure(1));
    writeVideo(MovAvi,F);
end
close(MovAvi)


stop

%%

% qv = load('Data2D/QvYHalf_99.dat');
% movies = size(qv,1)/(nx*m)
% qv = reshape(qv,nx,m,movies); %This are already in m/s
% 
% qr = load('Data2D/QrYHalf_99.dat');
% qr = reshape(qr,nx,m,movies); %This are already in m/s
% 
% qt = qv+qr;
% 
% figure(1)
% hold on
% for mov=1:movies
%     figure(1)
%     clf
%     contour(x,z,transpose(qt(:,:,mov)),100); colorbar;
%     print('Frame','-depsc',figure(1))
% end
% stop

%%


% qv3D =  load('Data3D/Qv_09.dat');
% qr3D =  load('Data3D/Qr_09.dat');
% 
% qt3D = qv3D+qr3D;
% 
% max(max(max(qv3D)))
% qt3D = reshape(qt3D,nx,ny,m);
% value = 0.1*max(max(max(qt3D)));
% isosurface(x,y,z,qt3D,value)
% axis([0 Lx 0 Ly 0 Lz])
% value
% stop

%%

% qv_ini = load('Data1D/Qv0_99.dat');
% qv_ini = reshape(qv_ini,m,1); %This are already in m/s
% 
% qvs = load('Data1D/Qvs_99.dat');
% qvs = reshape(qvs,m,1); %This are already in m/s
% 
% qvAnom = load('Data2D/QvYHalf_99.dat');
% movies = size(qvAnom,1)/(nx*m)
% qvAnom = reshape(qvAnom,nx,m,movies); %This are already in m/s
% 
% qr = load('Data2D/QrYHalf_99.dat');
% movies = size(qr,1)/(nx*m)
% qr = reshape(qr,nx,m,movies); %This are already in m/s
% 
% qvAnom_mean = zeros(m,movies);
% qr_mean = zeros(m,movies);
% for mov=1:movies
%     for iz = 1:m
%         qvAnom_mean(iz,mov) = max(qvAnom(:,iz,mov))/1;
%         qr_mean(iz,mov) = max(qr(:,iz,mov))/1;
%     end
% end
% 
% figure(1)
% hold on
% mov=movies;
% plot(qv_ini+qvAnom_mean(:,mov)+qr_mean(:,mov),z,'b -')
% plot(qvs,z,'r -')
% 
% stop
% 
% mov =1;
% max(max(abs(qvAnom(:,:,mov))))
% contour(x,z,transpose(qvAnom(:,:,mov))); colorbar;
% 
% stop

%%

% %Rain water:
% 
% qv_mov = load('Data3D/Qv_09.dat');
% size(qv_mov)/(nx*ny*m)
% qv_mov=reshape(qv_mov,nx,ny,m);
% 
% qr_mov = load('Data3D/Qr_09.dat');
% size(qr_mov)/(nx*ny*m)
% qr_mov=reshape(qr_mov,nx,ny,m);
% 
% qt_mov = qv_mov*0+qr_mov;
% 
% %Rain water 3D View
% xslice=Lx/2; %0:dx:Lx-dx; %10*(Lx-dx);% 10*Lx/2; %
% yslice=0:dy:(Ly-dy); %10*(Ly-dy); %0:10*Ly/32:10*(Ly-dy); %10*Lx/2; %
% zslice=Lz/2; %0:dz:Lz; 
% clf
% figure(1)
% %subplot(2,1,1)
% %h=slice(10*xGrid,10*yGrid,10*zGrid,(epsilon^2)*qr_mov,yslice,xslice,zslice);
% h=slice(y,x,z,qt_mov,yslice,xslice,zslice);
% xlabel('y','FontSize',15)
% ylabel('x','FontSize',15)
% set(gca,'XDir','rev')
% zlabel('z','FontSize',15)
% %title(['Squall line at T=', num2str(Tmov),' hours'],'FontSize',15)
% %caxis([.1 9])
% %axis([0 200 800 1000 0 15])
% %axesm('mercator','Origin',[0 128 0])
% %view(108,48)
% view(133,32)
% %axis tight %xy %equal %auto %manual
% alpha('color')
% set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
% c=colorbar;
% width1=get(gca,'position');
% width=get(c,'Position');
% width(3)=0.01;
% set(c,'Position',width)
% set(gca,'position',width1)
% %print('Qr3DT9','-depsc',figure(1))
% stop


%%

% qv_ini = load('Data1D/Qv0_99.dat');
% qv_ini = reshape(qv_ini,m,1); %This are already in m/s
% 
% qvs = load('Data1D/Qvs_99.dat');
% qvs = reshape(qvs,m,1); %This are already in m/s
% 
% 
% qv3D =  load('Data3D/Qv_01.dat');
% qv3D = reshape(qv3D,nx,ny,m);
% 
% qvZmax = zeros(m,1);
% for iz=1:m
%     qvZmax(iz,1) = sum(sum(qv3D(:,:,iz)))/nx/ny;
% end
% 
% hold on
% plot(qv_ini+qvZmax(:),z,'b -')
% plot(qvs,z,'r -')
% stop


%%

% qv0 = load('Data1D/Qv0_99.dat');
% qv0 = reshape(qv0,m,1); %This are already in m/s
% 
% qvs = load('Data1D/Qvs_99.dat');
% qvs = reshape(qvs,m,1); %This are already in m/s
% 
% figure(1)
% hold on
% plot(qv0,z,'b -',qvs,z,'r --')
% stop

%%

% movies=1183; %6*60+1;
% 
% qrZavg = load('Data2D/QrZavg_99.dat'); %load('ThetaZpt5Km_99.dat'); %
% moviesqr=round(size(qrZavg,1)/(nx*ny))-1;
% qrZavg2 = qrZavg(1:moviesqr*nx*ny,1);
% qrZavg2=reshape(qrZavg2,nx,ny,moviesqr);
% 
% % writerObj = VideoWriter('QrZavg.avi');
% % writerObj.FrameRate = 5;
% % open(writerObj);
% for mov=0*60+1:moviesqr
%     
%     Tmov=NoStepsMov*(mov-1);
% 
%     clf
%     figure(1)
%     contour(x,y,transpose(qrZavg2(:,:,mov)),30); colorbar;
%     caxis([0.1 3])
%     
%     axis([0 (Lx-dx) 0 (Ly-dy)])
%     axis on
%     xlabel('x in km');
%     ylabel('y in km');    
%     
%     title(['z-avg qr (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
%     
% %     frame = getframe(figure(1));
% %     writeVideo(writerObj,frame);
% end
% %close(writerObj);
% stop

%%

% movies=6*60+1;
% ubar = load('Ubar_99.dat'); 
% size(ubar)/(m)
% ubar=reshape(ubar,m,movies);
% 
% figure(1)
% hold on
% plot(ubar(:,1),z,'b -','LineWidth',2)
% plot(ubar(:,1*60+1),z,'b --','LineWidth',2)
% plot(ubar(:,2*60+1),z,'b .','LineWidth',2)
% plot(ubar(:,3*60+1),z,'black -','LineWidth',2)
% plot(ubar(:,4*60+1),z,'black --','LineWidth',2)
% plot(ubar(:,5*60+1),z,'black .','LineWidth',2)
% plot(ubar(:,6*60+1),z,'r --','LineWidth',2)
% xlabel('$\bar{u}$','interpreter','latex','FontSize',20)
% ylabel('z (km)','FontSize',20)
% h=legend('FARE: $\bar u$(T=0 hrs)','FARE: $\bar u$(T=1 hrs)','FARE: $\bar u$(T=2 hrs)','FARE: $\bar u$(T=3 hrs)','FARE: $\bar u$(T=4 hrs)','FARE: $\bar u$(T=5 hrs)','FARE: $\bar u$(T=6 hrs)');
% set(h,'interpreter','latex');
% set(h,'FontSize',10)
% print('UbarFARE3DT0to6','-depsc',figure(1))
% stop

%%


movies=6*60+1;
ubar = load('Vbar_99.dat'); 
size(ubar)/(m)
ubar=reshape(ubar,m,movies);

figure(1)
hold on
plot(ubar(:,1),z,'b -','LineWidth',2)
plot(ubar(:,1*60+1),z,'b --','LineWidth',2)
plot(ubar(:,2*60+1),z,'b .','LineWidth',2)
plot(ubar(:,3*60+1),z,'black -','LineWidth',2)
plot(ubar(:,4*60+1),z,'black --','LineWidth',2)
plot(ubar(:,5*60+1),z,'black .','LineWidth',2)
plot(ubar(:,6*60+1),z,'r --','LineWidth',2)
xlabel('$\bar{u}$','interpreter','latex','FontSize',20)
ylabel('z (km)','FontSize',20)
h=legend('FARE: $\bar v$(T=0 hrs)','FARE: $\bar v$(T=1 hrs)','FARE: $\bar v$(T=2 hrs)','FARE: $\bar v$(T=3 hrs)','FARE: $\bar v$(T=4 hrs)','FARE: $\bar v$(T=5 hrs)','FARE: $\bar v$(T=6 hrs)');
set(h,'interpreter','latex');
set(h,'FontSize',10)
set(h,'Location','NorthWest')
print('VbarFARE3DT0to6','-depsc',figure(1))
stop

%%

movies=6*60+1;
NoStepsMov=1/60;
vbar=load('vbar_99.dat');
size(vbar)/(m)
vbar=reshape(vbar,m,movies);

MovAvi1 = avifile('vbar.avi');
MovAvi1.fps = 10;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov; %In real time (hours)

    clf
    figure(1)
    hold on
    plot(vbar(1:m-1,1),z_u(1:m-1),'b -','LineWidth',2);
    plot(vbar(1:m-1,mov),z_u(1:m-1),'r .','LineWidth',2);
    axis([-1 5 0 15])
    xlabel('u bar');
    ylabel('z');    
    
    title(['T = ', num2str(Tmov,'% 10.2f'),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end

MovAvi1 = close(MovAvi1);
stop

%%

qrYHalf = load('qrYHalf_99.dat');
size(qrYHalf)/(nx*m)
qrYHalf=reshape(qrYHalf,nx,m,movies);

 
MovAvi1 = avifile('QrYHalf.avi');
MovAvi1.fps = 1;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov; %In real time (hours)

    clf
    figure(1)
    contour(x,z,transpose(qrYHalf(:,:,mov)),30); colorbar;
    caxis([0.3 15])
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qr (x,y) contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end
MovAvi1 = close(MovAvi1);

stop

%%

MovAvi1 = avifile('QvZMax.avi');
MovAvi1.fps = 1;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,(epsilon^2)*transpose(qvZMax(:,:,mov)),30); colorbar;
    %caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qv (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end

MovAvi1 = close(MovAvi1);


%%%%%%%%Water varpo Yhalf

MovAvi1 = avifile('QvYHalf.avi');
MovAvi1.fps = 1;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(qvYHalf(:,:,mov)),30); colorbar;
    caxis([0 5*10^(-3)])
    %caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qv (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end

MovAvi1 = close(MovAvi1);

stop 

%%%%%%%%%%theta at low levels

MovAvi1 = avifile('ThetaZpt5Km.avi');
MovAvi1.fps = 1;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,transpose(3*thetaZpt5Km(:,:,mov)),30); colorbar;
    caxis([-3 0.5])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['theta (x,y) contour at z=0.5 km in Kelvin, T = ', num2str(Tmov),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end

MovAvi1 = close(MovAvi1);

%%%%Steady squall line

MovAvi1 = avifile('QrYHalfSteady.avi');
MovAvi1.fps = 1;
for mov=it_0:it_1
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(qrYHalfSteady(:,:,mov)),30); colorbar;
    caxis([0.3*10^(-3) 10^(-3)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['Steady qr (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
    F = getframe(figure(1));
    MovAvi1 = addframe(MovAvi1,F);
end

MovAvi1 = close(MovAvi1);





