function pl = plotCircle3D(center,normal,radius,color,linewidth, grids)

theta=linspace(0,2*pi,grids);
v=null(normal');
points=repmat(center,1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),zeros(1,length(points(2,:))),points(3,:),color,'LineWidth',linewidth);
pl = plot3(points(1,:),points(2,:),points(3,:),color,'LineWidth',linewidth);

end