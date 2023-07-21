ps = Xs(1:3,:);
p = Xint(1:3,:);
pf = Xf(1:3,:);
qs = Xs(4:7,:);
q=Xint(4:7,:);
qf=Xf(4:7,:);
vs = Xs(8:10,:);
v = Xint(8:10,:);
vf = Xf(8:10,:);
ws = Xs(11:end,:);
w = Xint(11:end,:);
wf = Xf(11:end,:);

ind_convert = 1:13;
ind_Fc = ind_convert;
ind_Fc(admissible_control) = [];
ind_F = [ind_Fc admissible_control];

control_ind = zeros(1,13);
for i = 1:13
   control_ind(ind_F(i)) = i;
end

control = U(control_ind,:);

vd = control(8:10,:);
wd = control(11:13,:);

% plot cost vs s
figure('units','normalized','Position',[.25 .25 .5 .5])
plot(s(2:end),cost(2:end));
hold on
title('cost vs s')
plot(s(2:end),cost_all(2:end),'*m');
plot(s(2:end),cost_u(2:end),'r');
legend('total cost','cost of all inputs', 'cost of admissible inputs')

% position ,orientation, and corresponding controls
figure('units','normalized','outerposition',[0 0 1 1],'Name','Position and Orientation')
subplot(6,2,1)
plot(t,ps(1,:),'r', t,pf(1,:),'*g', t,p(1,:),'b');
ylabel('$x$','Interpreter','latex')
subplot(6,2,3)
plot(t,ps(2,:),'r', t,pf(2,:),'*g', t,p(2,:),'b');
ylabel('$y$','Interpreter','latex')
subplot(6,2,5)
plot(t,ps(3,:),'r', t,pf(3,:),'*g', t,p(3,:),'b');
ylabel('$z$','Interpreter','latex')
subplot(6,2,7)
plot(t,vs(1,:),'r', t,vf(1,:),'*g', t,v(1,:),'b');
ylabel('$v_x$','Interpreter','latex')
subplot(6,2,9)
plot(t,vs(2,:),'r', t,vf(2,:),'*g', t,v(2,:),'b');
ylabel('$v_x$','Interpreter','latex')
subplot(6,2,11)
plot(t,vs(3,:),'r', t,vf(3,:),'*g', t,v(3,:),'b');
ylabel('$v_x$','Interpreter','latex')


subplot(6,2,2)
plot(t,control(1,:));
ylabel('$v_x$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(6,2,4)
plot(t,control(2,:));
ylabel('$v_y$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(6,2,6)
plot(t,control(3,:));
ylabel('$v_z$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(6,2,8)
plot(t,vd(1,:));
ylabel('$\dot{v}_x$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(6,2,10)
plot(t,vd(2,:));
ylabel('$\dot{v}_y$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(6,2,12)
plot(t,vd(3,:));
ylabel('$\dot{v}_z$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );


% linear & angular velocity , and corresponding controls
figure('units','normalized','outerposition',[0 0 1 1],'Name','Linear and Augular Velocity')
subplot(7,2,1)
plot(t,qs(1,:),'r', t,qf(1,:),'*g', t,q(1,:),'b');
ylabel('$q_0$','Interpreter','latex')
subplot(7,2,3)
plot(t,qs(2,:),'r', t,qf(2,:),'*g', t,q(2,:),'b');
ylabel('$q_1$','Interpreter','latex')
subplot(7,2,5)
plot(t,qs(3,:),'r', t,qf(3,:),'*g', t,q(3,:),'b');
ylabel('$q_2$','Interpreter','latex')
subplot(7,2,7)
plot(t,qs(4,:),'r', t,qf(4,:),'*g', t,q(4,:),'b');
ylabel('$q_3$','Interpreter','latex')
subplot(7,2,9)
plot(t,ws(1,:),'r', t,wf(1,:),'*g', t,w(1,:),'b');
ylabel('$\omega_x$','Interpreter','latex')
subplot(7,2,11)
plot(t,ws(2,:),'r', t,wf(2,:),'*g', t,w(2,:),'b');
ylabel('$\omega_y$','Interpreter','latex')
subplot(7,2,13)
plot(t,ws(3,:),'r', t,wf(3,:),'*g', t,w(3,:),'b');
ylabel('$\omega_z$','Interpreter','latex')



subplot(7,2,2)
plot(t,control(4,:));
ylabel('$\dot{q}_0$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,4)
plot(t,control(5,:));
ylabel('$\dot{q}_1$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,6)
plot(t,control(6,:));
ylabel('$$\dot{q}_2$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,8)
plot(t,control(7,:));
ylabel('$$\dot{q}_3$ virtual','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,10)
plot(t,wd(1,:));
ylabel('$\dot{\omega}_x$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,12)
plot(t,wd(2,:));
ylabel('$\dot{\omega}_y$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
subplot(7,2,14)
plot(t,wd(3,:));
ylabel('$\dot{\omega}_z$','Interpreter','latex')
hold on
plot( [t(1) t(end)],0*ones(2,1) );
hold on
%%

qnorm = zeros(1,length(t));
qsnorm = zeros(1,length(t));
qfnorm = zeros(1,length(t));
for i = 1:length(t)
    qnorm(i) = norm(q(:,i));
    qsnorm(i) = norm(qs(:,i));
    qfnorm(i) = norm(qf(:,i));
end
figure('Name','Check unit quaternion','units','normalized','Position',[.25 .25 .5 .5])
plot(t,qsnorm,'r', t,qfnorm,'*g', t,qnorm,'b');
%% homotopy plot for norm
figure('units','normalized','Position',[.25 .25 .5 .5]);

qnorm_homo = zeros(1,length(t));
for i = 1:length(t)
    qnorm_homo(i) = norm( squeeze(sol(1,i,4:7)) );
end
plot(t,ones(1,length(t)), 'y', 'LineWidth',2);
%axis([0 t(end) min(qnorm_homo)-0.01 max(qnorm)+0.01])
xlabel('time, s')
ylabel('q norm')
hold on
grid on
plt = plot(t,qnorm_homo,'k', 'LineWidth',2);
plot(t,qnorm_homo, '-.r', 'LineWidth',2);
legend('norm = 1','q norm','norm of initial guess')

pause;

for j = 1:sgrids
    for i = 1:length(t)
        qnorm_homo(i) = norm( squeeze(sol(j,i,4:7)) );
    end
    plt.YDataSource='qnorm_homo';
    refreshdata(plt,'caller');
    
    pause(0.02);
end

%% x y z homotopy animation

figure('units','normalized','Position',[.25 .25 .5 .5]);

x=sol(1,:,1);
y=sol(1,:,2);
z=sol(1,:,3);
h1=plot3(x,y,z,'k','LineWidth',2);
axis([min(p(1,:))-0.2 max(p(1,:))+0.2 min(p(2,:))-0.2 max(p(2,:))+0.2 min(p(3,:))-0.2 max(p(3,:))+0.2
    ]);
grid ON;
title('x-y-z configuration space curve');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
pause;

for i=1:sgrids
    x=sol(i,:,1);h1.XDataSource='x';
    y=sol(i,:,2);h1.YDataSource='y';
    z=sol(i,:,3);h1.ZDataSource='z';
    
    refreshdata(h1,'caller');
    drawnow;
end
%% animation
figure('units','normalized','Position',[.25 .25 .5 .5]);

R = quat2rotm(q');
Rs = quat2rotm(qs');
vlength = 0.2;
D_box = [0.5 0.3 0.05];

R_axis1 = [p(:,1) p(:,1)+vlength*R(:,1,1)];
R_axis2 = [p(:,1) p(:,1)+vlength*R(:,2,1)];
R_axis3 = [p(:,1) p(:,1)+vlength*R(:,3,1)];
R_plot = plot3(R_axis1(1,:),R_axis1(2,:),R_axis1(3,:),'r',...
               R_axis2(1,:),R_axis2(2,:),R_axis2(3,:),'g',...
               R_axis3(1,:),R_axis3(2,:),R_axis3(3,:),'b');
R_plot(1).LineWidth = 2;
R_plot(2).LineWidth = 2;
R_plot(3).LineWidth = 2;
box on
grid on
hold on
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis([min(p(1,:))-D_box(1)/2-0.1 max(p(1,:))+D_box(1)/2+0.1 ...
      min(p(2,:))-D_box(2)/2-0.1 max(p(2,:))+D_box(2)/2+0.1 ...
      min(p(3,:))-D_box(3)/2-0.1 max(p(3,:))+D_box(3)/2+0.1]);


Rs_axis1 = [ps(:,1) ps(:,1)+vlength*Rs(:,1,1)];
Rs_axis2 = [ps(:,1) ps(:,1)+vlength*Rs(:,2,1)];
Rs_axis3 = [ps(:,1) ps(:,1)+vlength*Rs(:,3,1)];
Rs_plot = plot3(Rs_axis1(1,:),Rs_axis1(2,:),Rs_axis1(3,:),':c',...
               Rs_axis2(1,:),Rs_axis2(2,:),Rs_axis2(3,:),':c',...
               Rs_axis3(1,:),Rs_axis3(2,:),Rs_axis3(3,:),':c');
Rs_plot(1).LineWidth = 2;
Rs_plot(2).LineWidth = 2;
Rs_plot(3).LineWidth = 2;

xyz_plot = plot3(p(1,:),p(2,:),p(3,:),'k','LineWidth',2);
plot3(ps(1,:),ps(2,:),ps(3,:),':c','LineWidth',2);

legend([R_plot; xyz_plot; Rs_plot(1)],{'ax','ay','az','integrated xyz path','AGHF solution'})
%legend([R_plot; xyz_plot; Rs_plot(1)],{'ax','ay','az','integrated xyz path','AGHF solution'},'AutoUpdate','off')

box_plt = DrawBox(p(:,1),R(:,:,1),D_box,'m',0.3);


pause;

for i = 1:length(t)
    
    delete(box_plt);
    
    R_axis1 = [p(:,i) p(:,i)+vlength*R(:,1,i)];
    R_axis2 = [p(:,i) p(:,i)+vlength*R(:,2,i)];
    R_axis3 = [p(:,i) p(:,i)+vlength*R(:,3,i)];
    
    Rs_axis1 = [ps(:,i) ps(:,i)+vlength*Rs(:,1,i)];
    Rs_axis2 = [ps(:,i) ps(:,i)+vlength*Rs(:,2,i)];
    Rs_axis3 = [ps(:,i) ps(:,i)+vlength*Rs(:,3,i)];
    
    R_plot(1).XDataSource='R_axis1(1,:)';R_plot(1).YDataSource='R_axis1(2,:)';R_plot(1).ZDataSource='R_axis1(3,:)';
    R_plot(2).XDataSource='R_axis2(1,:)';R_plot(2).YDataSource='R_axis2(2,:)';R_plot(2).ZDataSource='R_axis2(3,:)';
    R_plot(3).XDataSource='R_axis3(1,:)';R_plot(3).YDataSource='R_axis3(2,:)';R_plot(3).ZDataSource='R_axis3(3,:)';
    Rs_plot(1).XDataSource='Rs_axis1(1,:)';Rs_plot(1).YDataSource='Rs_axis1(2,:)';Rs_plot(1).ZDataSource='Rs_axis1(3,:)';
    Rs_plot(2).XDataSource='Rs_axis2(1,:)';Rs_plot(2).YDataSource='Rs_axis2(2,:)';Rs_plot(2).ZDataSource='Rs_axis2(3,:)';
    Rs_plot(3).XDataSource='Rs_axis3(1,:)';Rs_plot(3).YDataSource='Rs_axis3(2,:)';Rs_plot(3).ZDataSource='Rs_axis3(3,:)';
    
    refreshdata(R_plot,'caller');
    refreshdata(Rs_plot,'caller');
    
    box_plt = DrawBox(p(:,i),R(:,:,i),D_box,'m',0.3);
    
    drawnow;
    
    t_pause = 0.1/length(t);
    pause(t_pause);
end


