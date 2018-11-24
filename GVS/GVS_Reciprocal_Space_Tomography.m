%% Read tomography data

me=load('Cylinder_mesh_dense_cart.mat','X_cyl_dense','Y_cyl_dense','Z_cyl_dense');
Int=load('GVS_volumetric_cyl_sm.mat','I_cyl_cylint');

Xc=me.X_cyl_dense/100; %convert to nm^-1
Yc=me.Y_cyl_dense/100;
Zc=me.Z_cyl_dense/100;

I=Int.I_cyl_cylint;

%Rotate grid to compensate offset in polar angle
[Phic,Rc,Zc]=cart2pol(Xc,Yc,Zc);
Phic=mod(Phic*180/pi-80,360)*pi/180; %100° rotation
[Xc,Yc,Zc]=pol2cart(Phic,Rc,Zc);

Xc=permute(Xc,[2,1,3]);
Yc=permute(Yc,[2,1,3]);
Zc=permute(Zc,[2,1,3]);
I=permute(I,[2,1,3]);


%% Transform to Cartesian grid
x=-0.5:0.01:0.5;
y=x; z=x;
[Xg,Yg,Zg]=meshgrid(x,y,z);


%I_cart=griddata(Xc,Yc,Zc,I,Xg,Yg,Zg);
%save('GVS_Cart_intensity_m80.mat','I_cart')
Int2=load('GVS_Cart_intensity_m80.mat','I_cart');
I_cart_GVS=Int2.I_cart;
I_cart_GVS=permute(I_cart_GVS,[2,1,3]);

Rmax=0.3; %forward scattered intensity
I_cart_GVS(Xg.^2+Yg.^2+Zg.^2<Rmax^2)=NaN; %masking noise within the sphere with r<0.03


%% - Transform coordinate system from [11-2,111,1-10]to [100,010,001]
T=[1/sqrt(6)*[1;1;-2],-1/sqrt(3)*[1;1;1],1/sqrt(2)*[1;-1;0]];
A=[Xg(:),Yg(:),Zg(:)]*T';
Xt=reshape(A(:,1),size(Xg));
Yt=reshape(A(:,2),size(Yg));
Zt=reshape(A(:,3),size(Zg));
%save('Cartesian_mesh_Base.mat','Xt','Yt','Zt')

figure(12)
wx=9.67;
wy=9.22;
thr=2.5;
cla
hold on
%Plot_Slicing_Planes;
Plot_RST(Xt,Yt,Zt,I_cart_GVS,thr);
%daspect([1 1 1]) 
view(100,17);
camup([1,-1,0])
lim=0.5;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square

%export_fig('GVS_RCT_3D','-png','-m4')

figure(13)
Pos=[3.7888    2.7093    9.7155    8.0010];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_cart_GVS,thr);
%view(1/sqrt(3)*[-1,-1,1])

set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
lim=0.5;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
camup([1,-1,0])
axis square
axis tight
grid off

%%Rotation
alpha=0:0.1:6.28;
grid off
axis on
for i=1:numel(alpha)
    vw=cos(alpha(i))*[0,0,1]+1/sqrt(2)*sin(alpha(i))*[1,1,0];
    view(vw+[0,0,0])
   
    axis vis3d
    camup([1,-1,0])
    drawnow
    %pause(0.1)
    export_fig(['RST_rotation_base',num2str(i),'.png'],'-png','-nocrop')
end

%
   
figure(14)
Pos=[13.5467    2.6670    9.1863    8.0645];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_cart_GVS,thr);
%view(-54.7,0);
view([0,0,1])
camup([1,-1,0])
lim=0.5;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square
grid off
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);



figure(15)
Pos=[22.7753    2.7093    6.8792    8.0645];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_cart_GVS,thr);
%view(+35.2,0);
view([1,1,0])
camup([1,-1,0])
lim=0.5;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square
grid off
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
        


figure(16)
Pos=[29.6968    2.0320   10.3717    9.3133];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_cart_GVS,thr);
view([1,1,-2])
camup([1,-1,0])
lim=0.5;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
%axis square
grid off

%{
figure(33)
thr_low=3;
thr_high=7;

ind=find(I<thr_high & I>thr_low);
[K,L,M]=ind2sub([361,61,128],ind);

vx=[];
vy=[];
vz=[];
for i=1:length(K)
    vx(i)=Xc(K(i),L(i),M(i));
    vy(i)=Yc(K(i),L(i),M(i));
    vz(i)=Zc(K(i),L(i),M(i));
end

subplot(2,3,1)
cla
box on
grid off
scatter3(vx,vy,vz,1,I(ind))
phi=54.7*pi/180;
view(0,0);

subplot(2,3,2)
cla
box on
grid off
scatter3(vx,vy,vz,1,I(ind))
view(54.7,0);


subplot(2,3,3)
cla
box on
grid off
scatter3(vx,vy,vz,1,I(ind))
view(-35.2,0);
hsp2 = get(gca, 'Position');  
%{
%% Plot cross sections
colormap(jet)
res=1;
offset=-9.5;

subplot(2,3,4)
cla
phi=(offset-70.52)*pi/180;
z1=-50; z2=50;
r1=-50; r2=50;
[X2D,Z2D,I_100]= Slicing_plane_intensity(Xc,Yc,Zc,I,z1,z2,r1,r2,res,phi);
pcolor(X2D,Z2D,I_100); shading flat;
axis square

subplot(2,3,5)
cla
phi=(offset-35.2)*pi/180;
z1=-50; z2=50;
r1=-50; r2=50;
[X2D,Z2D,I_110]= Slicing_plane_intensity(Xc,Yc,Zc,I,z1,z2,r1,r2,res,phi);
pcolor(X2D,Z2D,I_110); shading flat;
axis square

subplot(2,3,6)
phi=(offset+54.7)*pi/180;
z1=-50; z2=50;
r1=-50; r2=50;
[X2D,Z2D,I_100]= Slicing_plane_intensity(Xc,Yc,Zc,I,z1,z2,r1,r2,res,phi);
pcolor(X2D,Z2D,I_100); shading flat;
axis square
%}



figure(34)
cla

scatter3(vx,vy,vz,2,I(ind))
daspect([0.75 0.75 1])
%view([1,1,1])
hold on

%% Plot slicing planes
offset=0;
phi=(offset+54.7)*pi/180;
[Xs,Ys,Zs]=Slicing_plane(-0.5,0.5,0.5,phi);
surf(squeeze(Xs),squeeze(Ys),squeeze(Zs),mean(I(ind)));

phi=(offset-70.52)*pi/180;
[Xs,Ys,Zs]=Slicing_plane(-0.5,0.5,0.5,phi);
surf(squeeze(Xs),squeeze(Ys),squeeze(Zs),mean(I(ind)));

phi=(offset-35.2)*pi/180;
[Xs,Ys,Zs]=Slicing_plane(-0.5,0.5,0.5,phi);
surf(squeeze(Xs),squeeze(Ys),squeeze(Zs),mean(I(ind)));
alpha 0.2





%rotating in 3D
%{
alpha=0:0.1:6.28;
grid off
axis on
for i=1:numel(alpha)
    vw=cos(alpha(i))*[1,0,0]+sin(alpha(i))*[0,1,0];
    view(vw+[0,0,0])
    axis vis3d
    %drawnow
    %pause(0.1)
    %export_fig(['RST_rotation_front',num2str(i),'.png'],'-png','-nocrop')
end
%}




%{
hold on
%v1=[1,0,0];
%v2=[0,1,0];
%v3=[0,0,1];
n=[1,1,1];

xp=-70:1:70;
yp=-70:1:70;
[Xpl,Ypl]=meshgrid(xp,yp);
H=0;
Zpl=(H-n(1)*Xpl-n(2)*Ypl)/n(3);
%surf(Xpl,Ypl,Zpl);
caxis([thr_low,thr_high])

%Ipl=griddata(Xc,Yc,Zc,I,Xpl,Ypl,Zpl);

Xq=reshape(Xpl,1,[]);
Yq=reshape(Ypl,1,[]);
Zq=reshape(Zpl,1,[]);
%Ipl=interp3(Xc,Yc,Zc,I,Xq,Yq,Zq);
%Ipl=reshape(Ipl,size(Zpl));
%}
%{
hslice = surf(linspace(-rmax,rmax,100),...
   linspace(zmin,zmax,100),...
   zeros(100));

rotate(hslice,[-1,0,0],-45)

caxis([thr_low,thr_high])
xd = get(hslice,'XData');
yd = get(hslice,'YData');
zd = get(hslice,'ZData');
%}
%}