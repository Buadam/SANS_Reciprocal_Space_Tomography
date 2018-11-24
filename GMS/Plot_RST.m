function []= Plot_RST(Xc,Yc,Zc,I,thr,col,alpha)
%Coordinate system
%p = patch(isosurface(Xc,Yc,Zc,floor(I_cart/1.5),1));

if nargin<7
    alpha=1;
end
if nargin<6
    col=[0.4,0.4,1];
end

I2=zeros(size(I));

%thr=3;
I2(I>thr)=100;
reshape(I2,size(I));

p = patch(isosurface(Xc,Yc,Zc,I2,99));
%isonormals(Xc,Yc,Zc,I2,p)
p.FaceColor = col;
p.FaceAlpha=alpha;
p.EdgeColor = [0.3,0.3,0.3];
p.EdgeColor = 'none';%[0.3,0.3,0.3];
p.LineWidth=0.5;
p.EdgeAlpha=0.6;

p.FaceLighting = 'flat';
p.AmbientStrength = 0.3;
p.DiffuseStrength = 0.4;
p.SpecularStrength = 1;
p.SpecularExponent = 8;
p.BackFaceLighting = 'reverselit';

light('Position',[0.8 0 0.8],'Style','local')
lighting flat
grid on
box on
%xlabel('q_x (nm^{-1})')
%ylabel('q_y (nm^{-1})')
%zlabel('q_z (nm^{-1})')

set(gcf, 'Color', 'w');
set(gca,'FontSize',12);

end