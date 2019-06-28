
function [Q_111,Q_1m1m1,Q_m11m1,Q_m1m11]= Model_Modulated_rings(alpha,qmax)

if nargin<2
    qmax=1;
end
if nargin<1
    alpha=0.2;
end

phi=linspace(0,2*pi,1000)';
%theta=pi/2+alpha*sqrt(2)/3*sin(3*phi);
theta=pi/2+sqrt(2)/3*alpha/(1+alpha)*sin(3*phi);
qx=qmax*sin(theta).*cos(phi);
qy=qmax*sin(theta).*sin(phi);
qz=qmax*cos(theta);
Q=[qx,qy,qz];

A1=[1/sqrt(2)*[1;-1;0],1/sqrt(6)*[1;1;-2],1/sqrt(3)*[1;1;1]];
A2=[1/sqrt(2)*[0;-1;1],-1/sqrt(6)*[2;1;1],1/sqrt(3)*[1;-1;-1]];
A3=[1/sqrt(2)*[1;0;-1],-1/sqrt(6)*[1;2;1],1/sqrt(3)*[-1;1;-1]];
A4=[-1/sqrt(2)*[1;-1;0],-1/sqrt(6)*[1;1;2],1/sqrt(3)*[-1;-1;1]];


Q_111=Q*A1';
Q_1m1m1=Q*A2';
Q_m11m1=Q*A3';
Q_m1m11=Q*A4';

%{
figure(4324)
cla
hold all
grid on

plot3(Q_111(:,1),Q_111(:,2),Q_111(:,3))
plot3(Q_1m1m1(:,1),Q_1m1m1(:,2),Q_1m1m1(:,3))
plot3(Q_m11m1(:,1),Q_m11m1(:,2),Q_m11m1(:,3))
plot3(Q_m1m11(:,1),Q_m1m11(:,2),Q_m1m11(:,3))
axis square
axis vis3d
%}
end
%Q=[1;1;1]*cos(theta)*sqrt(3)+sin(theta)*([1;-1;0]*cos(phi)*1/sqrt(2)+[1;1;-2]*1/sqrt(6));

