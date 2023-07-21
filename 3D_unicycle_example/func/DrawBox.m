function plt = DrawBox(c,R,D,color,FaceAlpha)
% given center and rotation matrix, dimensions, draw box
p1 = c + R*[-D(1) -D(2) -D(3)]'/2;
p2 = c + R*[D(1) -D(2) -D(3)]'/2;
p3 = c + R*[D(1) D(2) -D(3)]'/2;
p4 = c + R*[-D(1) D(2) -D(3)]'/2;
p5 = c + R*[-D(1) -D(2) D(3)]'/2;
p6 = c + R*[D(1) -D(2) D(3)]'/2;
p7 = c + R*[D(1) D(2) D(3)]'/2;
p8 = c + R*[-D(1) D(2) D(3)]'/2;

f1 = [p1 p2 p3 p4];
f2 = [p5 p6 p7 p8];
f3 = [p1 p2 p6 p5];
f4 = [p3 p4 p8 p7];
f5 = [p2 p3 p7 p6];
f6 = [p1 p4 p8 p5];

X = [f1(1,:)' f2(1,:)' f3(1,:)' f4(1,:)' f5(1,:)' f6(1,:)'];
Y = [f1(2,:)' f2(2,:)' f3(2,:)' f4(2,:)' f5(2,:)' f6(2,:)'];
Z = [f1(3,:)' f2(3,:)' f3(3,:)' f4(3,:)' f5(3,:)' f6(3,:)'];
plt = fill3(X,Y,Z,color,'FaceAlpha',FaceAlpha);

end