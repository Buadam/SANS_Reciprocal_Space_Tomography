omega_deg = sym(linspace(0,360,361));
omega_rad = sym(omega_deg*pi/180);
alpha=-0.2;

% Lab -> Experimental coordinate trafo
A_1m10 = sym([1/sqrt(2), 0, 1/sqrt(2); 1/sqrt(2), 0, -1/sqrt(2); 0,1,0]);

% Lab -> Rhombohedral coordinate trafo
A_1m1m1 = sym([1/sqrt(2), 1/sqrt(6), 1/sqrt(3); 1/sqrt(2), -1/sqrt(6), -1/sqrt(3); 0,2/sqrt(6), -1/sqrt(3)]);

R= (A_1m1m1)'*A_1m10;
% x axis rotated by omega in the rh coordinate system
r_omega_rh = R*[cos(omega_rad); sin(omega_rad); zeros(size(omega_rad))];   % matrix of size (3, omega)

% parametrize phi(omega) -- offset by 90° so that parametrization starts
% from the orthogonal direction to the x axis (omega=0)
phi = atan2(r_omega_rh(2,:), r_omega_rh(1,:))+pi/2;  % parameter vector of size (omega)

% Calculate parametric curve in the rh coord system
theta = pi/2 + sqrt(2)/3*alpha/(alpha+1)*sin(3*phi);

% Calculate the q-vector space curve in Descartes coordinates and transform
% it to the experimental coordinate system
Q = R'*[sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];
chi_rad = acos(Q(3,:));
% chi_rad = acos(Q(3,:)./(sqrt(Q(1,:).^2 + Q(2,:).^2 + Q(3,:).^2)));
% chi_rad = atan((sqrt(Q(1,:).^2 + Q(2,:).^2))./ Q(3,:));
chi_deg = chi_rad * 180/pi;

figure(1)
subplot(1,3,1)
plot(omega_deg-35.2, phi*180/pi) 
subplot(1,3,2)
plot(omega_deg-35.2, theta*180/pi)
subplot(1,3,3)
plot(omega_deg-35.2, chi_deg) 

figure(2)
cla
X=sin(chi_rad).*cos(omega_rad);
Y=sin(chi_rad).*sin(omega_rad);
Z=cos(chi_rad);
scatter3(Q(1,:), Q(2,:), Q(3,:))
hold on
scatter3(X, Y, Z)
axis square
box on

