function [x,z,I_0_sm]=SANS_image_conditioning(I_0,c_x,c_y,pix_x,pix_y,Rmax,sigma,moving)

%% Define meshgrid with the center in the origin
ratio=pix_x/pix_y;
x=0:1:pix_x-1; x=x-c_x; %coordinate of the last pixel = pix_x-1
z=0:1:pix_y-1; z=ratio*(z-c_y);
[X,Z]=meshgrid(x,z);


%% Convolutional Smoothing
if sigma~=0
    sz = 2*ceil(2.6 * sigma) + 1; 
    mask = fspecial('gauss', sz, sigma);
end
I_0_sm=[];

%% Smooth all layers and cut forward scattered intensity with NaN-s
rg=X.^2+Z.^2; %circle around the center with rg radius square
for i=1:size(I_0,3)
    if sigma ~=0
        I_temp = conv2(I_0(:,:,i), mask, 'same');
    else
        I_temp=I_0(:,:,i);
    end
    I_temp(rg<Rmax^2)=NaN; %cut forward scattering
    I_0_sm(:,:,i)=I_temp;
end

%% Moving window to average consecutive layers
if moving~=0
    I_0_sm=movmean(I_0_sm,moving,3); %central average for the phi angles. the size does not change
end


%{
%% Crop data to [-xmax:xmax], [-ymax,ymax]
if (x_max ~=0 && y_max~=0)
    ind1=find(z>-y_max,1);
    ind2=find(z>y_max,1);
    ind3=find(x>-x_max,1);
    ind4=find(x>x_max,1);

    x2=x(ind3:ind4-1);
    z2=z(ind1:ind2-1);
    [X2,Z2]=meshgrid(x2,z2);
    I_0_crop=I_0_sm(ind1:ind2-1, ind3:ind4-1,:);
end
%}