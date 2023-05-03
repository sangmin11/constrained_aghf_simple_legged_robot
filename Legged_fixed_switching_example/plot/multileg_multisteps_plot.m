xs=Xs(1,:);
ys=Xs(2,:);
ths=Xs(3,:);
xds=Xs(4,:);
yds=Xs(5,:);
thds=Xs(6,:);

x=Xint(1,:);
y=Xint(2,:);
th=Xint(3,:);
xd=Xint(4,:);
yd=Xint(5,:);
thd=Xint(6,:);


xf=Xf(1,:);
yf=Xf(2,:);
thf=Xf(3,:);
xdf=Xf(4,:);
ydf=Xf(5,:);
thdf=Xf(6,:);

v1=U(1,:);
v2=U(2,:);
v3=U(3,:);
v4=U(4,:);
v5=U(5,:);
v6=U(6,:);


Force = [];
Pos = [];
dPos = [];
Force_s = [];
Pos_s = [];
Force_f = [];
Pos_f = [];
u_leg = [];

for i=1:k_leg
    ind = 7+4*(i-1);
    Force = [Force; Xint(ind:ind+1,:)];
    Pos = [Pos; Xint(ind+2:ind+3,:)];
    Force_s = [Force_s; Xs(ind:ind+1,:)];
    Pos_s = [Pos_s; Xs(ind+2:ind+3,:)];
    Force_f = [Force_f; Xf(ind:ind+1,:)];
    Pos_f = [Pos_f; Xf(ind+2:ind+3,:)];
    u_leg = [u_leg; U(ind:ind+3,:)];
    dPos = [dPos; U(ind+2:ind+3,:)];
end


% plot cost vs s
figure
plot(s(2:end),cost(2:end));
hold on
title('cost vs s')
plot(s(2:end),cost_all(2:end),'*m');
plot(s(2:end),cost_u(2:end),'r');
legend('total cost','cost of all inputs', 'cost of admissible inputs')

figure('units','normalized','outerposition',[0 0 1 1])

subplot(6,2,2)
plot(t,v1,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_1$','Interpreter','latex')
subplot(6,2,4)
plot(t,v2,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_2$','Interpreter','latex')
subplot(6,2,6)
plot(t,v3,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_3$','Interpreter','latex')
subplot(6,2,8)
plot(t,v4,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_4$','Interpreter','latex')
subplot(6,2,10)
plot(t,v5,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_5$','Interpreter','latex')
subplot(6,2,12)
plot(t,v6,[t(1) t(end)],0*ones(2,1),'r');
ylabel('$v_6$','Interpreter','latex')

subplot(6,2,1)
plot(t,xs,'r', t,xf,'*g', t,x,'b');
ylabel('$x$','Interpreter','latex')
subplot(6,2,3)
plot(t,ys,'r', t,yf,'*g', t,y,'b');
ylabel('$y$','Interpreter','latex')
subplot(6,2,5)
plot(t,ths,'r', t,thf,'*g', t,th,'b');
ylabel('$\theta$','Interpreter','latex')
subplot(6,2,7)
plot(t,xds,'r', t,xdf,'*g', t,xd,'b');
ylabel('$\dot{x}$','Interpreter','latex')
subplot(6,2,9)
plot(t,yds,'r', t,ydf,'*g', t,yd,'b');
ylabel('$\dot{y}$','Interpreter','latex')
subplot(6,2,11)
plot(t,thds,'r', t,thdf,'*g', t,thd,'b');
ylabel('$\dot{\theta}$','Interpreter','latex')

figure('units','normalized','outerposition',[0 0 1 1])

for i=1:k_leg
    
    ind = 2*(i-1)+1;
    ind4 = 4*(i-1)+1;
    ind8 = 8*(i-1)+1;
    
    subplot(4*k_leg,2,ind8)
    plot(t,Force_s(ind,:),'r', t,Force_f(ind,:),'*g', t,Force(ind,:),'b');
    ylabel(['$F' num2str(i) '_x$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+2)
    plot(t,Force_s(ind+1,:),'r', t,Force_f(ind+1,:),'*g', t,Force(ind+1,:),'b');
    ylabel(['$F' num2str(i) '_y$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+4)
    plot(t,Pos_s(ind,:),'r', t,Pos_f(ind,:),'*g', t,Pos(ind,:),'b');
    ylabel(['$p' num2str(i) '_x$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+6)
    plot(t,Pos_s(ind+1,:),'r', t,Pos_f(ind+1,:),'*g', t,Pos(ind+1,:),'b');
    ylabel(['$p' num2str(i) '_y$'],'Interpreter','latex')
    
    subplot(4*k_leg,2,ind8+1)
    plot(t,u_leg(ind4,:),[t(1) t(end)],0*ones(2,1),'r');
    ylabel(['$\dot{F}' num2str(i) '_x$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+3)
    plot(t,u_leg(ind4+1,:),[t(1) t(end)],0*ones(2,1),'r');
    ylabel(['$\dot{F}' num2str(i) '_y$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+5)
    plot(t,u_leg(ind4+2,:),[t(1) t(end)],0*ones(2,1),'r');
    ylabel(['$\dot{p}' num2str(i) '_x$'],'Interpreter','latex')
    subplot(4*k_leg,2,ind8+7)
    plot(t,u_leg(ind4+3,:),[t(1) t(end)],0*ones(2,1),'r');
    ylabel(['$\dot{p}' num2str(i) '_y$'],'Interpreter','latex')
end

%% animation
figure
plot(xs,ys,'LineWidth',3);
hold on;
plot(x,y,'m','LineWidth',2)
xterrain = linspace(min(min(Pos(1:2:2*k_leg-1,:))),max(max(Pos(1:2:2*k_leg-1,:))),100);
yterrain = terrain_func(xterrain);
plot(xterrain,yterrain,'-.k','LineWidth',2)
%force_scale = 0.001;
%h_force=plot(Xforce, Yforce,'b','LineWidth',2);
h_force=quiver(Pos(1:2:2*k_leg-1,1), Pos(2:2:2*k_leg,1),force_scale*Force(1:2:2*k_leg-1,1), force_scale*Force(2:2:2*k_leg,1),0, 'linewidth',2,'Color','b', 'MaxHeadSize',0.5);
legend('HF solution','integrated path','terrain','contact force')
scale = 0.2;
X1=[cos(ths); sin(ths)];
X2=[-sin(ths); cos(ths)];
X11=[cos(th); sin(th)];
X22=[-sin(th); cos(th)];
Xplot1 = [x(1) x(1)+scale*X1(1,1)];
Yplot1 = [y(1) y(1)+scale*X1(2,1)];
Xplot2 = [x(1) x(1)+scale*X2(1,1)];
Yplot2 = [y(1) y(1)+scale*X2(2,1)];
Xplot11 = [xdrift(1,1) xdrift(1,1)+scale*X11(1,1)];
Yplot11 = [xdrift(2,1) xdrift(2,1)+scale*X11(2,1)];
Xplot22 = [xdrift(1,1) xdrift(1,1)+scale*X22(1,1)];
Yplot22 = [xdrift(2,1) xdrift(2,1)+scale*X22(2,1)];

h1=plot(Xplot1, Yplot1,'r','LineWidth',3);
h2=plot(Xplot2, Yplot2,'g','LineWidth',3);
h11=plot(Xplot11, Yplot11,'--r','LineWidth',2);
h22=plot(Xplot22, Yplot22,'--g','LineWidth',2);

box_c = [xdrift(1,1) xdrift(2,1)]';
Rot = [cos(th(1)) -sin(th(1)); sin(th(1)) cos(th(1))];
box1 = box_c + Rot*[w_robot/2 h_robot/2]';
box2 = box_c +Rot*[-w_robot/2 h_robot/2]';
box3 = box_c +Rot*[-w_robot/2 -h_robot/2]';
box4 = box_c +Rot*[w_robot/2 -h_robot/2]';
BOX = [box1 box2 box3 box4 box1];
box = plot(BOX(1,:),BOX(2,:),'c','LineWidth',2);

plegx = [];
plegy = [];
for i=1:k_leg
    pknee = find_knee([x(1);y(1)],Pos( (2*(i-1)+1):(2*(i-1)+2) ,1),l1,l2,1);
    plegx = [plegx [x(1);pknee(1);Pos(2*(i-1)+1 ,1)] ];
    plegy = [plegy [y(1);pknee(2);Pos(2*(i-1)+2 ,1)] ];
end

legs_plot = plot(plegx,plegy,'c','LineWidth',2);

axis equal
%pbaspect([1 1 1]);
axis([min(xs)-0.4, max(xs)+0.4, -0.4, max(ys)+0.4]);
grid ON;
xlabel('x');
ylabel('y');
pause;


for i = 1:length(t)
   
    Xplot1=[xs(i) xs(i)+scale*X1(1,i)];h1.XDataSource='Xplot1';
    Yplot1=[ys(i) ys(i)+scale*X1(2,i)];h1.YDataSource='Yplot1';
    Xplot2=[xs(i) xs(i)+scale*X2(1,i)];h2.XDataSource='Xplot2';
    Yplot2=[ys(i) ys(i)+scale*X2(2,i)];h2.YDataSource='Yplot2';
    Xplot11=[x(i) x(i)+scale*X11(1,i)];h11.XDataSource='Xplot11';
    Yplot11=[y(i) y(i)+scale*X11(2,i)];h11.YDataSource='Yplot11';
    Xplot22=[x(i) x(i)+scale*X22(1,i)];h22.XDataSource='Xplot22';
    Yplot22=[y(i) y(i)+scale*X22(2,i)];h22.YDataSource='Yplot22';
    set(h_force,'xdata',Pos(1:2:2*k_leg-1,i),'ydata',Pos(2:2:2*k_leg,i),'udata',force_scale*Force(1:2:2*k_leg-1,i),'vdata',force_scale*Force(2:2:2*k_leg,i));
    box_c = [x(i) y(i)]';
    Rot = [cos(th(i)) -sin(th(i)); sin(th(i)) cos(th(i))];
    box1 = box_c + Rot*[w_robot/2 h_robot/2]';
    box2 = box_c +Rot*[-w_robot/2 h_robot/2]';
    box3 = box_c +Rot*[-w_robot/2 -h_robot/2]';
    box4 = box_c +Rot*[w_robot/2 -h_robot/2]';
    BOX = [box1 box2 box3 box4 box1];
    box.XDataSource='BOX(1,:)';
    box.YDataSource='BOX(2,:)';
    
    plegx = [];
    plegy = [];
    for ii=1:k_leg
        pknee = find_knee([x(i);y(i)],Pos( (2*(ii-1)+1):(2*(ii-1)+2) ,i),l1,l2,1);
        plegx = [plegx [x(i);pknee(1);Pos(2*(ii-1)+1 ,i)] ];
        plegy = [plegy [y(i);pknee(2);Pos(2*(ii-1)+2 ,i)] ];
        legs_plot(ii).XDataSource='plegx(:,ii)';
        legs_plot(ii).YDataSource='plegy(:,ii)';
        refreshdata(legs_plot(ii),'caller');
    end
    
    
    %refreshdata;
    
    refreshdata(h1,'caller');
    refreshdata(h2,'caller');
    refreshdata(h11,'caller');
    refreshdata(h22,'caller');
    refreshdata(h_force,'caller');
    refreshdata(box,'caller');
    %refreshdata(legs_plot,'caller');
    
    drawnow;

    
    t_pause = 2/length(t);
    pause(t_pause);
   
end

