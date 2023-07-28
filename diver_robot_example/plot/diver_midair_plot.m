ths=Xs(1,:);
q1s=Xs(2,:);
q2s=Xs(3,:);
thds=Xs(4,:);
q1ds=Xs(5,:);
q2ds=Xs(6,:);

th=Xint(1,:);
q1=Xint(2,:);
q2=Xint(3,:);
thd=Xint(4,:);
q1d=Xint(5,:);
q2d=Xint(6,:);

thf=Xf(1,:);
q1f=Xf(2,:);
q2f=Xf(3,:);
thdf=Xf(4,:);
q1df=Xf(5,:);
q2df=Xf(6,:);

u1=U(1,:);
u2=U(2,:);
u3=U(3,:);
u4=U(4,:);
u5=U(5,:);
u6=U(6,:);

bu1=bufull(1,:);
bu2=bufull(2,:);
bu3=bufull(3,:);
bu4=bufull(4,:);
bu5=bufull(5,:);
bu6=bufull(6,:);


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
plot(t,u1,tint,bu1);
hold on
subplot(6,2,4)
plot(t,u2,tint,bu2);
hold on
subplot(6,2,6)
plot(t,u3,tint,bu3);
hold on
subplot(6,2,8)
plot(t,u4,tint,bu4);
hold on
subplot(6,2,10)
plot(t,u5,tint,bu5);
hold on
subplot(6,2,12)
plot(t,u6,tint,bu6);
hold on


subplot(6,2,1)
plot(t,ths,'r', t,thf,'*g', t,th,'b');
ylabel('$\theta$','Interpreter','latex')
legend('AGHF solution','integrated path','full path(with inadmissible control)')
subplot(6,2,3)
plot(t,q1s,'r', t,q1f,'*g', t,q1,'b');
ylabel('$q_1$','Interpreter','latex')
subplot(6,2,5)
plot(t,q2s,'r', t,q2f,'*g', t,q2,'b');
ylabel('$q_2$','Interpreter','latex')
subplot(6,2,7)
plot(t,thds,'r', t,thdf,'*g', t,thd,'b');
ylabel('$\dot{\theta}$','Interpreter','latex')
subplot(6,2,9)
plot(t,q1ds,'r', t,q1df,'*g', t,q1d,'b');
ylabel('$\dot{q}_1$','Interpreter','latex')
subplot(6,2,11)
plot(t,q2ds,'r', t,q2df,'*g', t,q2d,'b');
ylabel('$\dot{q}_2$','Interpreter','latex')


%% animation
figure

vxc0 = -6;
vyc0 = g*(T/2);
xc0 = 0;
yc0 = 1.5;
Xcom = xc0 + vxc0*t;
Ycom = yc0 +(vyc0+vyc0-g*t)/2.*t; 
Pc = [Xcom;Ycom];
plot(Xcom,Ycom,'r')
%legend('CoM')
xlabel('x (m)');
ylabel('y (m)');
hold on

x = zeros(length(t),1);
y = zeros(length(t),1);
x1_b = -l1/2*sin(q1f(1));
y1_b =  l/2 + l1/2*cos(q1f(1));
x2_b =  l2/2*sin(q2f(1));
y2_b =  -l/2 - l2/2*cos(q2f(1));
R = [cos(thf(1))      -sin(thf(1));
     sin(thf(1))       cos(thf(1))];
X1 = [x(1) y(1)]' + R*[x1_b y1_b]';
X2 = [x(1) y(1)]' + R*[x2_b y2_b]'; 

Xc = ( m*[x(1) y(1)]' + m1*X1 + m2*X2 )/ (m+m1+m2);

pb1 = [x(1) y(1)]' + R*[0 l/2]' - Xc + Pc(:,1);
pb2 = [x(1) y(1)]' + R*[0 -l/2]' - Xc + Pc(:,1);
link_base=plot([pb1(1) pb2(1)], [pb1(2) pb2(2)],'c','LineWidth',2);
grid on
axis equal
%axis([-2 2 -2 2]);
axis([min(Xcom)-1 max(Xcom)+1 min(Ycom)-1.5 max(Ycom)+1.5]);
hold on


p1x_b = -l1*sin(q1f(1));
p1y_b =  l/2 + l1*cos(q1f(1));
p2x_b =  l2*sin(q2f(1));
p2y_b =  -l/2 - l2*cos(q2f(1));
p1 = [x(1) y(1)]' + R*[p1x_b p1y_b]' - Xc + Pc(:,1);
p2 = [x(1) y(1)]' + R*[p2x_b p2y_b]' - Xc + Pc(:,1);
p1x = [pb1(1) p1(1)];
p1y = [pb1(2) p1(2)];
p2x = [pb2(1) p2(1)];
p2y = [pb2(2) p2(2)];
link1=plot(p1x, p1y,'b','LineWidth',2);
link2=plot(p2x, p2y,'g','LineWidth',2);
CoM = plot(Pc(1,1),Pc(2,1),'ro','MarkerFaceColor','r','MarkerSize',8);
pause;

N_start = 20;
for i = 1:N_start
    frame(i) = getframe(gcf);
end
for i = 1:length(t)
    R = [cos(thf(i))      -sin(thf(i));
             sin(thf(i))       cos(thf(i))];
    x1_b = -l1/2*sin(q1f(i));
    y1_b =  l/2 + l1/2*cos(q1f(i));
    x2_b =  l2/2*sin(q2f(i));
    y2_b =  -l/2 - l2/2*cos(q2f(i));
    X1 = [x(i) y(i)]' + R*[x1_b y1_b]';
    X2 = [x(i) y(i)]' + R*[x2_b y2_b]';
    Xc = ( m*[x(i) y(i)]' + m1*X1 + m2*X2 )/ (m+m1+m2);
    pb1 = [x(i) y(i)]' + R*[0 l/2]' -Xc + Pc(:,i);
    pb2 = [x(i) y(i)]' + R*[0 -l/2]' -Xc + Pc(:,i);
    pbx = [pb1(1) pb2(1)]';
    pby = [pb1(2) pb2(2)]';
    link_base.XDataSource='pbx';
    link_base.YDataSource='pby';
    
    p1x_b = -l1*sin(q1f(i));
    p1y_b =  l/2 + l1*cos(q1f(i));
    p2x_b =  l2*sin(q2f(i));
    p2y_b =  -l/2 - l2*cos(q2f(i));
    p1 = [x(i) y(i)]' + R*[p1x_b p1y_b]' - Xc  + Pc(:,i);
    p2 = [x(i) y(i)]' + R*[p2x_b p2y_b]' - Xc  + Pc(:,i);
    p1x = [pb1(1) p1(1)];
    p1y = [pb1(2) p1(2)];
    p2x = [pb2(1) p2(1)];
    p2y = [pb2(2) p2(2)];
    link1.XDataSource='p1x';
    link1.YDataSource='p1y';
    link2.XDataSource='p2x';
    link2.YDataSource='p2y';
    
    xcom = Pc(1,i);
    ycom = Pc(2,i);
    CoM.XDataSource='xcom';
    CoM.YDataSource='ycom';
    
    refreshdata(link_base,'caller');
    refreshdata(link1,'caller');
    refreshdata(link2,'caller');
    refreshdata(CoM,'caller');
    
    drawnow;
    
    if( (i==1) || (mod(i,20)==0) || (i==length(t)) )
        plot([pb1(1) pb2(1)], [pb1(2) pb2(2)],'--c','LineWidth',2);
        plot(p1x, p1y,'--b','LineWidth',2);
        plot(p2x, p2y,'--g','LineWidth',2);
        plot(Pc(1,i),Pc(2,i),'ro','MarkerFaceColor','r','MarkerSize',8);
    end
    
    frame(i+N_start) = getframe(gcf);
    
    %t_pause = 5/length(t);
    %pause(t_pause);
end
N_end = 20;
for i = 1:N_end
    frame(i+N_start+length(t)) = getframe(gcf);
end

% video = VideoWriter('current_video.avi');
% open(video);

% writeVideo(video,frame);
% close(video);