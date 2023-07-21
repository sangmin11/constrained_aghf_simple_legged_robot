ps = Xs(1:3,:);
p = Xint(1:3,:);
pf = Xf(1:3,:);
qs = Xs(4:7,:);
q=Xint(4:7,:);
qf=Xf(4:7,:);


[N,Mt] = size(Xint);
ind_convert = 1:N;
ind_Fc = ind_convert;
ind_Fc(admissible_control) = [];
ind_F = [ind_Fc admissible_control];

control_ind = zeros(1,N);
for i = 1:N
   control_ind(ind_F(i)) = i;
end

control = U(control_ind,:);

%vd = control(8:10,:);
%wd = control(11:13,:);

%w0 = U(4,:);
%w = U([5 6 7],:);
%v = U([1 2 3],:);
w0 = control(4,:);
w = control(5:7,:);
v = control(1:3,:);

% plot cost vs s
figure
plot(s(2:end),cost(2:end));
hold on
title('cost vs s')
plot(s(2:end),cost_all(2:end),'*m');
plot(s(2:end),cost_u(2:end),'r');
legend('total cost','cost of all inputs', 'cost of admissible inputs')


figure('units','normalized','outerposition',[0 0 1 1],'Name','states and inputs')
subplot(7,2,1)
plot(t,ps(1,:),'r', t,pf(1,:),'*g', t,p(1,:),'b');
ylabel('$x$','Interpreter','latex')
subplot(7,2,3)
plot(t,ps(2,:),'r', t,pf(2,:),'*g', t,p(2,:),'b');
ylabel('$y$','Interpreter','latex')
subplot(7,2,5)
plot(t,ps(3,:),'r', t,pf(3,:),'*g', t,p(3,:),'b');
ylabel('$z$','Interpreter','latex')
subplot(7,2,7)
plot(t,qs(1,:),'r', t,qf(1,:),'*g', t,q(1,:),'b');
ylabel('$q_0$','Interpreter','latex')
subplot(7,2,9)
plot(t,qs(2,:),'r', t,qf(2,:),'*g', t,q(2,:),'b');
ylabel('$q_1$','Interpreter','latex')
subplot(7,2,11)
plot(t,qs(3,:),'r', t,qf(3,:),'*g', t,q(3,:),'b');
ylabel('$q_2$','Interpreter','latex')
subplot(7,2,13)
plot(t,qs(4,:),'r', t,qf(4,:),'*g', t,q(4,:),'b');
ylabel('$q_3$','Interpreter','latex')


subplot(7,2,2)
plot(t,v(1,:));
ylabel('$v_x$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,4)
plot(t,v(2,:));
ylabel('$v_y$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,6)
plot(t,v(3,:));
ylabel('$v_z$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,8)
plot(t,w0);
ylabel('$w_0$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,10)
plot(t,w(1,:));
ylabel('$\omega_x$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,12)
plot(t,w(2,:));
ylabel('$\omega_y$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,14)
plot(t,w(3,:));
ylabel('$\omega_z$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
hold on

%% x y z homotopy animation

figure;

x=sol(1,:,1);
y=sol(1,:,2);
z=sol(1,:,3);
h1=plot3(x,y,z,'k','LineWidth',2);
axis equal
axis([-1 1 -1 1 -1 1]);
%axis([min(p(1,:))-0.1 max(p(1,:))+0.1 min(p(2,:))-0.1 max(p(2,:))+0.1 min(p(3,:))-0.1 max(p(3,:))+0.1]);
grid ON;
title('x-y-z configuration space curve');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
pause;

%{
N_start = 30;
for i = 1:N_start
    frame(i) = getframe(gcf);
end
%}
for i=1:sgrids
    x=sol(i,:,1);h1.XDataSource='x';
    y=sol(i,:,2);h1.YDataSource='y';
    z=sol(i,:,3);h1.ZDataSource='z';
    
    refreshdata(h1,'caller');
    drawnow;
    
    %frame(i+N_start) = getframe(gcf);
    
    t_pause = 0.2/length(t);
    pause(t_pause);
    
end

%{
video = VideoWriter('current_video.avi');
video.FrameRate = round(length(t)/2);
video.Quality = 50;
open(video);

writeVideo(video,frame);
close(video);
%}
%% animation
figure('units','normalized','outerposition',[0 0 1 1]);

R = quat2rotm(q');
Rs = quat2rotm(qs');

Rcircle = 0.01;
linewidth = 2;
circle_grids = 20;
path_plot = plotCircle3D( p(:,1), R(:,:,1)*[1 0 0]', Rcircle, 'b', linewidth, circle_grids );
hold on
box on
grid on
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis([min(p(1,:))-0.1 max(p(1,:))+0.1 ...
      min(p(2,:))-0.1 max(p(2,:))+0.1 ...
      min(p(3,:))-0.1 max(p(3,:))+0.1]);
paths_plot = plotCircle3D( ps(:,1), Rs(:,:,1)*[1 0 0]', Rcircle, ':c', linewidth, circle_grids );

d_base = 0.05;
patch([0 0 0 0],[-d_base d_base d_base -d_base],[d_base d_base -d_base -d_base],'y');

p_top = zeros(3,length(t));
ps_top = zeros(3,length(t));
for i = 1:length(t)
    
    plotCircle3D( p(:,i), R(:,:,i)*[1 0 0]', Rcircle, 'b', linewidth, circle_grids );
    plotCircle3D( ps(:,i), Rs(:,:,i)*[1 0 0]', Rcircle, ':c', linewidth, circle_grids );
    p_top(:,i) = p(:,i)+R(:,:,i)*(Rcircle*[0 0 1]');
    ps_top(:,i) = ps(:,i)+Rs(:,:,i)*(Rcircle*[0 0 1]');
    
end

plot3(p_top(1,:),p_top(2,:),p_top(3,:),'g');
plot3(ps_top(1,:),ps_top(2,:),ps_top(3,:),':m');

patch([min(xlim) min(xlim) max(xlim) max(xlim)],[0 0 0 0],[min(zlim) max(zlim) max(zlim) min(zlim)],'r','FaceAlpha',.1,'EdgeAlpha',0);
%patch('Faces',Fpatch,'Vertices',Vpatch,'FaceVertexCData',Cpatch,'FaceColor','flat','FaceAlpha',.2,'EdgeAlpha',0);
if( ~isempty(fixed_points) )
    for i = 1:length(fixed_points(1,:))
        plot3([fixed_points(1,i) fixed_points(1,i)],[min(ylim) max(ylim)],[fixed_points(3,i) fixed_points(3,i)],'r','LineWidth',2);
    end
    plot3(p_top(1,:),p_top(2,:),p_top(3,:),'g');
    marker_plot = scatter3(fixed_points(1,:),fixed_points(2,:),fixed_points(3,:),100,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
    legend([path_plot; paths_plot; marker_plot],{'integrated path','AGHF solution','Markers'})
else
    legend([path_plot; paths_plot],{'integrated path','AGHF solution'})
end

