
load('Cartesian_mesh_Lab.mat') %load standard base mesh
%load('Cartesian_mesh_Base.mat') %load standard base mesh
%load('GMS_intensity_raw.mat','I_cart_cartint') %load raw data
%load('GMS_intensity_gauss_31_3.mat','I_cart') %load smoothed data
load('GMS_intensity_fullsymmetrized.mat','I_av') %load fully symmetrized data

I_fit=I_av; %select data to fit

I_fit(isnan(I_fit))=0;

thr_low=4;
ind=find(I_fit>thr_low);
[K,L,M]=ind2sub([171,171,171],ind);
Qx=squeeze(Xc(1,L,1)); %ndgrid: 1st index K, meshgrid: 2nd index K
Qy=squeeze(Yc(K,1,1)); %ndgrid: 2nd index L, meshgrid: 1st index L
Qz=squeeze(Zc(1,1,M));

%% Transform coordinate system from [11-2,111,1-10] to [100,010,001]
%%
T=[1/sqrt(6)*[1;1;-2],1/sqrt(3)*[1;1;1],1/sqrt(2)*[1;-1;0]];
A=[Xc(:),Yc(:),Zc(:)]*T';
Xt=reshape(A(:,1),size(Xc));
Yt=reshape(A(:,2),size(Yc));
Zt=reshape(A(:,3),size(Zc));

Q=[Qx',Qy,Qz];

Qexp=Q*T';
W=I_fit(ind); %weights of the experimental q-vectors, i.e. SANS counts at the given vectors

%% Visualize Q vectors that have a larger intensity than thr_low
%{
figure(43242)
cla
grid on
hold all
scatter3(Qexp(:,1),Qexp(:,2),Qexp(:,3),2,I_fit(ind),'MarkerFaceAlpha',.4)
axis square
%}

%% Fitting rings to data
%{
alpha0=-0.22;
qmax0=0.65;
p0=[alpha0,qmax0];
options = optimoptions(@fminunc,'Display','iter');
modelfun=@(p)Rings_error(p,Qexp,I_fit,ind); %model function is weighted likewise
[params,fval,exitflag,output,grad,hessian]=fminunc(modelfun, p0, options);

alpha_Fit=params(1) %fit:-0.156
qmax_Fit=params(2) %fit: 0.644
err = sqrt(diag(inv(hessian)))
%}
%% Define the four rhombohedral axes
%==========================================================================
n1=1/sqrt(3)*[1;1;1];
n2=1/sqrt(3)*[-1;-1;1];
n3=1/sqrt(3)*[-1;1;-1];
n4=1/sqrt(3)*[1;-1;-1];


%% Separate rhombohedral domains
Qp=[(Qexp*n1).^2,(Qexp*n2).^2,(Qexp*n3).^2,(Qexp*n4).^2]; 
[~,ind2]=min(Qp,[],2); %find Q vectors that minimize (ni*Q)^2 for i=1,2,3,4

Q_111=[Qexp(ind2==1,1),Qexp(ind2==1,2),Qexp(ind2==1,3)];
w_111=W(ind2==1);
Q_m1m11=[Qexp(ind2==2,1),Qexp(ind2==2,2),Qexp(ind2==2,3)];
w_m1m11=W(ind2==2);
Q_m11m1=[Qexp(ind2==3,1),Qexp(ind2==3,2),Qexp(ind2==3,3)];
w_m11m1=W(ind2==3);
Q_1m1m1=[Qexp(ind2==4,1),Qexp(ind2==4,2),Qexp(ind2==4,3)];
w_1m1m1=W(ind2==4);


%Plot the contributions of the four rhombohedral domains

figure(3633)
cla
grid on
hold all
r1=scatter3(Qexp(ind2==1,1),Qexp(ind2==1,2),Qexp(ind2==1,3),2,'red');
r2=scatter3(Qexp(ind2==2,1),Qexp(ind2==2,2),Qexp(ind2==2,3),2,'blue');
r3=scatter3(Qexp(ind2==3,1),Qexp(ind2==3,2),Qexp(ind2==3,3),2,'green');
r4=scatter3(Qexp(ind2==4,1),Qexp(ind2==4,2),Qexp(ind2==4,3),2,'yellow');
r1.MarkerEdgeAlpha=0.3;
r2.MarkerEdgeAlpha=0.3;
r3.MarkerEdgeAlpha=0.3;
r4.MarkerEdgeAlpha=0.3;

view(100,13);
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square




%Fitted parameters:
%alpha_Fit=-0.156: symmetrzied thr=4
%           -0.224: raw data, thr=3.6
%            -0.286: raw data, thr=1
%            -0.3: raw data, thr=0
%            -0.139: Gaussian, thr=0.5
%            -0.107: Gaussian, thr=1.7
%qmax_Fit=0.644: symmetrized thr=4
%           0.643: raw data, thr=3.6
%           0.68: raw data, thr=1
%           0.7: raw data, thr=0
%           0.7: Gaussian, thr=0.5
%           0.684: Gaussian, thr=1.7


alpha_Fit=-0.135;
qmax_Fit=0.643;

%Visualize Fit
[Qmod_111,Qmod_1m1m1,Qmod_m11m1,Qmod_m1m11]= Model_Modulated_rings(alpha_Fit,qmax_Fit);

figure(321312)

thr=3.6; %Intensity threshold for visualization : best option: 3.6
wx=9.67;
wy=9.22;

cla
hold on
%%Plot_Slicing_Planes;
Plot_RST(Xt,Yt,Zt,I_fit,thr);
set(gcf,'Units','centimeters');
%daspect([1 1 1]) 

plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')

set(gcf, 'Position', [1.1007 1.7357 wx wy]);
%view(-162,47); %view1
%view(117,52); %top view
view(100,13);
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square

beta=0:0.1:6.28;
grid off
axis on
for i=1:1%numel(beta)
    vw=cos(beta(i))*[1,1,1]/sqrt(3)+sin(beta(i))*[1,1,-2]/sqrt(6);
    view(vw+[0,0,0])
    camup([1,-1,0])
    
    axis vis3d
    drawnow
    %pause(0.1)
    %export_fig(['RST_rotation_model',num2str(i),'.png'],'-pdf','-nocrop')
end


%%Plot views
col=[0.4,0.4,1]; %blue
alpha=0.3;
thr=3.6; %Intensity threshold for visualization : best option: 3.6
figure(52)
wx=9.9;
wy=9.22;

cla
hold on
Plot_RST(Xt,Yt,Zt,I_av,thr,col,alpha);
plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')
set(gcf,'Units','centimeters');
set(gcf, 'Position', [1.1007 1.7357 wx wy]);
view(100,13);
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
axis square
grid on
export_fig('Figures/GMS_persp_fit.png','-pdf','-nocrop');



figure(53)
Pos=[5.9690    3.1750    8.6    8.4];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_av,thr,col,alpha);
plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')
view(1/sqrt(3)*[1,1,1])
camup([1,-1,0])
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis square
grid on
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
export_fig('Figures/GMS_111_fit.png','-pdf');


figure(54)
Pos=[14.1817    3.1327    8.6    8.4];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_av,thr,col,alpha);
plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')
view([0,0,1])
camup([1,-1,0])
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis square
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
%axis square
grid on
export_fig('Figures/GMS_001_fit.pdf');

figure(55)
Pos=[22.7753    2.7093    7.56    8.4];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_av,thr,col,alpha);
plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')
view(1/sqrt(2)*[1,1,0])
camup([1,-1,0])
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis square
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
%axis square
grid on
export_fig('Figures/GMS_110_fit.png','-pdf');



figure(56)
Pos=[29.9297    3.5348    8.6    8.4];
cla
hold on
Plot_RST(Xt,Yt,Zt,I_av,thr,col,alpha);
plot3(Qmod_111(:,1),Qmod_111(:,2),Qmod_111(:,3),'linewidth',2.5,'color','r')
plot3(Qmod_1m1m1(:,1),Qmod_1m1m1(:,2),Qmod_1m1m1(:,3),'linewidth',2.5,'color','y')
plot3(Qmod_m11m1(:,1),Qmod_m11m1(:,2),Qmod_m11m1(:,3),'linewidth',2.5,'color','b')
plot3(Qmod_m1m11(:,1),Qmod_m1m11(:,2),Qmod_m1m11(:,3),'linewidth',2.5,'color','g')
view(1/sqrt(6)*[1,1,2])
camup([1,-1,0])
lim=0.65;
xlim([-lim,lim])
ylim([-lim,lim])
zlim([-lim,lim])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis square
set(gcf,'Units','centimeters');
set(gcf, 'Position', Pos);
%axis square
grid on
export_fig('Figures/GMS_112_fit.png','-pdf');
%==========================================================================



%}
%{
%% ------------------Extract intensity along the rings-------------------%%
q1=0.4;
q2=0.9;
dq=0.05;
q_dist=q1:dq:q2;

%Convert mesh to vector components
thr_low=1;
ind=find(I_fit>thr_low);
[K,L,M]=ind2sub([171,171,171],ind);
Qx=squeeze(Xc(1,L,1)); %ndgrid: 1st index K, meshgrid: 2nd index K
Qy=squeeze(Yc(K,1,1)); %ndgrid: 2nd index L, meshgrid: 1st index L
Qz=squeeze(Zc(1,1,M));

%transform to [11-2;111;1-10] coord system
T=[1/sqrt(6)*[1;1;-2],1/sqrt(3)*[1;1;1],1/sqrt(2)*[1;-1;0]];
Q=[Qx',Qy,Qz];
Qall=Q*T';


I_ring_mean=zeros(size(q_dist));
I_ring_tot=zeros(size(q_dist));

for i=1:length(q_dist)
    [I_ring_mean,I_ring_tot(i),Q_ring{i}]=Get_intensity([alpha_Fit,q_dist(i)],Qall,dq,I_fit);
end

figure(43242)
cla
grid on
hold all
for i=1:length(q_dist)
    scatter3(Q_ring{i}(:,1),Q_ring{i}(:,2),Q_ring{i}(:,3),1,I_ring_mean(i)*ones(size(Q_ring{i},1),1),'MarkerFaceAlpha',.4)
end
axis square


figure(42342)
cla
box on
hold on
plot(q_dist,I_ring_mean./q_dist)


figure(42343)
cla
box on
hold on
plot(q_dist,I_ring_tot./q_dist)
%}