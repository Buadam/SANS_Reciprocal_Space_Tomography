function I_0_sm=SANS_image_conditioning(I_0,c_x,c_y,pix_x,pix_y,Rmax,sigma)


x=0.5-c_x:1:pix_x-c_x-0.5;
z=0.5-c_y:1:pix_y-c_y-0.5;
[X,Z]=meshgrid(x,z);

rg=X.^2+Z.^2;


sz = 2*ceil(2.6 * sigma) + 1; % See note below
mask = fspecial('gauss', sz, sigma);

I_0_sm=[];

for i=1:size(I_0,3)
    I_temp = conv2(I_0(:,:,i), mask, 'same');
    I_temp(rg<Rmax^2)=NaN;
    I_0_sm(:,:,i)=I_temp;
end

%{
figure(12)
cla
pcolor(X,Z,I_0_sm(:,:,4)); shading flat;
caxis([0,70])
axis square

hold on
plot(0.5,0.5,'rx')
%}
end