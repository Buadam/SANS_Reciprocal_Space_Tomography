%function [Q_111,Q_1m1m1,Q_m11m1,Q_m1m11]= Model_Modulated_rings_numerical(alpha)



%if nargin<1
%    alpha=-0.2;
%end
alpha=0.2;
n1=1/sqrt(3)*[1;1;1];
n2=1/sqrt(3)*[-1;-1;1];
n3=1/sqrt(3)*[-1;1;-1];
n4=1/sqrt(3)*[1;-1;-1];

%Create 2D grid on the surface of a sphere (|Q|=const)
r=[1,1.001];
phi=0:0.01:2*pi;
theta=0:0.01:pi;
[phi_g,theta_g,r]=meshgrid(phi,theta,r);

Qx=sin(theta_g).*cos(phi_g);
Qy=sin(theta_g).*sin(phi_g);
Qz=cos(theta_g);


U_111=(Qx*n1(1)+Qy*n1(2)+Qz*n1(3)).^2+abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_m1m11=(Qx*n2(1)+Qy*n2(2)+Qz*n2(3)).^2+abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_m11m1=(Qx*n3(1)+Qy*n3(2)+Qz*n3(3)).^2+abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_1m1m1=(Qx*n4(1)+Qy*n4(2)+Qz*n4(3)).^2+abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);

minU=min(min(min(U_111)));
eps=0.05;
ind1=find(U_111<minU+abs(eps*minU)); [K,L,M]=ind2sub(size(Qx),ind1);
Qnum_111=zeros(numel(ind1),3);
for i=1:numel(ind1)
    Qnum_111(i,:)=[Qx(K(i),L(i),M(i)),Qy(K(i),L(i),M(i)),Qz(K(i),L(i),M(i))];
end
Qnum_111

U_111_min=(Qx*n1(1)+Qy*n1(2)+Qz*n1(3)).^2-abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_m1m11_min=(Qx*n2(1)+Qy*n2(2)+Qz*n2(3)).^2-abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_m11m1_min=(Qx*n3(1)+Qy*n3(2)+Qz*n3(3)).^2-abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);
U_1m1m1_min=(Qx*n4(1)+Qy*n4(2)+Qz*n4(3)).^2-abs(alpha)*(Qx.^4+Qy.^4+Qz.^4);

minU=min(min(min(U_111_min)));
eps=0.05;
ind1=find(U_111_min<minU+abs(eps*minU)); [K,L,M]=ind2sub(size(Qx),ind1);
Qnum_111_min=zeros(numel(ind1),3);
for i=1:numel(ind1)
    Qnum_111_min(i,:)=[Qx(K(i),L(i),M(i)),Qy(K(i),L(i),M(i)),Qz(K(i),L(i),M(i))];
end
Qnum_111_min

%ind2=find(U_m1m11<eps*min(min(min(U_m1m11))));
%ind3=find(U_m11m1<eps*min(min(min(U_m11m1))));
%ind4=find(U_1m1m1<eps*min(min(min(U_1m1m1))));

%Q_1m1m1=Q*A2';
%Q_m11m1=Q*A3';
%Q_m1m11=Q*A4';

[Q_111,Q_1m1m1,Q_m11m1,Q_m1m11]= Model_Modulated_rings(abs(alpha),1);
[Q_111_min,Q_1m1m1_min,Q_m11m1_min,Q_m1m11_min]= Model_Modulated_rings(-abs(alpha),1);
[Q_111_0,Q_1m1m1_0,Q_m11m1_0,Q_m1m11_0]= Model_Modulated_rings(0,1);

%Value of the potential along the rings
U_Q=@(Q,n) (Q(:,1)*n(1)+Q(:,2)*n(2)+Q(:,3)*n(3)).^2+abs(alpha)*(Q(:,1).^4+Q(:,2).^4+Q(:,3).^4);
U_Q_min=@(Q,n) (Q(:,1)*n(1)+Q(:,2)*n(2)+Q(:,3)*n(3)).^2-abs(alpha)*(Q(:,1).^4+Q(:,2).^4+Q(:,3).^4);

UQ111=U_Q(Q_111,n1);
UQm111=U_Q_min(Q_111_min,n1);

phi=linspace(0,2*pi,1000)*180/pi;
p111=exp(-UQ111)./trapz(phi,exp(-UQ111)); p111=(p111-min(p111))/mean(p111);

%(Q_111(:,1)*n1(1)+Q_111(:,2)*n1(2)+Q_111(:,3)*n1(3)).^2+abs(alpha)*(Q_111(:,1).^4+Q_111(:,2).^4+Q_111(:,3).^4)
figure(4324)
wx=12;
wy=8;
set(gcf,'Units','centimeters');
set(gcf, 'Position', [1.12 11.34 wx wy]);

cla
hold all
grid on
scatter3(Qnum_111(:,1),Qnum_111(:,2),Qnum_111(:,3),3,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter3(Q_111(:,1),Q_111(:,2),Q_111(:,3),10,UQ111,'filled')
scatter3(Q_111_0(:,1),Q_111_0(:,2),Q_111_0(:,3),1,'k.')

axis square
view(65,10)
colormap(jet)

%{
p1 = patch(isosurface(Qx,Qy,Qz,floor(U_111*70),min(min(min(floor(U_111*70))))));
p1.EdgeColor='r';
p2 = patch(isosurface(Qx,Qy,Qz,floor(U_m1m11),min(min(min(U_m1m11)))));
p2.EdgeColor='b';
p3 = patch(isosurface(Qx,Qy,Qz,floor(U_m11m1),min(min(min(U_m11m1)))));
p3.EdgeColor='g';
p4 = patch(isosurface(Qx,Qy,Qz,floor(U_1m1m1),min(min(min(U_1m1m1)))));
p4.EdgeColor='y';

plot3(Q_111(:,1),Q_111(:,2),Q_111(:,3))
plot3(Q_1m1m1(:,1),Q_1m1m1(:,2),Q_1m1m1(:,3))
plot3(Q_m11m1(:,1),Q_m11m1(:,2),Q_m11m1(:,3))
plot3(Q_m1m11(:,1),Q_m1m11(:,2),Q_m1m11(:,3))

plot3(Q_1m1m1(:,1),Q_1m1m1(:,2),Q_1m1m1(:,3))
plot3(Q_m11m1(:,1),Q_m11m1(:,2),Q_m11m1(:,3))
plot3(Q_m1m11(:,1),Q_m1m11(:,2),Q_m1m11(:,3))

%}
%end

figure(4325)
wx=12;
wy=8;
set(gcf,'Units','centimeters');
set(gcf, 'Position', [13.12 11.34 wx wy]);

cla
hold all
grid on
scatter3(Qnum_111_min(:,1),Qnum_111_min(:,2),Qnum_111_min(:,3),3,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter3(Q_111_min(:,1),Q_111_min(:,2),Q_111_min(:,3),10,exp(-(UQm111-mean(UQm111))/(abs(mean(UQm111)))),'filled')
scatter3(Q_111_0(:,1),Q_111_0(:,2),Q_111_0(:,3),1,'k.')

axis square
view(65,10)
colormap(jet)


figure(4326)
wx=12;
wy=8;
set(gcf,'Units','centimeters');
set(gcf, 'Position', [25.12 11.34 wx wy]);

cla
hold all
grid on

scatter3(Qnum_111(:,1),Qnum_111(:,2),Qnum_111(:,3),3,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter3(Q_111(:,1),Q_111(:,2),Q_111(:,3),10,U_Q(Q_111,n1),'filled')

scatter3(Qnum_111_min(:,1),Qnum_111_min(:,2),Qnum_111_min(:,3),3,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter3(Q_111_min(:,1),Q_111_min(:,2),Q_111_min(:,3),10,U_Q_min(Q_111_min,n1),'filled')

axis square
view(65,10)
colormap(jet)

figure(53453)
cla
hold all
grid on
phi=linspace(0,2*pi,1000)*180/pi;
U_111_pert_plusalpha=U_Q(Q_111,n1);%(Q_111(:,1)*n1(1)+Q_111(:,2)*n1(2)+Q_111(:,3)*n1(3)).^2+abs(alpha)*(Q_111(:,1).^4+Q_111(:,2).^4+Q_111(:,3).^4);
U_111_pert_minusalpha=U_Q_min(Q_111_min,n1);%(Q_111_min(:,1)*n1(1)+Q_111_min(:,2)*n1(2)+Q_111_min(:,3)*n1(3)).^2-abs(alpha)*(Q_111_min(:,1).^4+Q_111_min(:,2).^4+Q_111_min(:,3).^4);
plot(phi,U_111_pert_plusalpha)
plot(phi,U_111_pert_minusalpha)
legend(['alpha= ' num2str(abs(alpha))], ['alpha= ' abs(num2str(-abs(alpha)))])

