alpha=-0.135;
dn=2; %[1-1-1] or [-1,1,-1]
%% 
%==================Define trasformation matrices=========================%
% Lab -> Experimental coordinate trafo
%A_1m10 = [1/sqrt(2), 0, 1/sqrt(2); 1/sqrt(2), 0, -1/sqrt(2); 0,1,0];
A_1m10 = [1/sqrt(3), 1/sqrt(6), 1/sqrt(2); 1/sqrt(3), 1/sqrt(6), -1/sqrt(2); 1/sqrt(3), -2/sqrt(6),0];

% Lab -> Rhombohedral coordinate trafo
%A_1m1m1 = [1/sqrt(2), 1/sqrt(6), 1/sqrt(3); 1/sqrt(2), -1/sqrt(6), -1/sqrt(3); 0,2/sqrt(6), -1/sqrt(3)];
if dn == 1
    A_rh = [1/sqrt(2), 1/sqrt(6), 1/sqrt(3); 1/sqrt(2), -1/sqrt(6), -1/sqrt(3); 0,2/sqrt(6), -1/sqrt(3)];
else
    A_rh = [1/sqrt(2), -1/sqrt(6), -1/sqrt(3); 0, -2/sqrt(6), 1/sqrt(3); -1/sqrt(2), -1/sqrt(6), -1/sqrt(3)];
end

%Rhombohedral -> Experimental coordinate trafo
R= (A_rh)'*A_1m10;

%%
% ============== Calculate space curve ============================%

%equally spaced parametrization in the polar angle
phi_deg = linspace(0,359,360);  %starting from 0deg <110> direction
phi_rad = phi_deg*pi/180;

% Calculate parametric curve in the rh coord system
theta = pi/2 + sqrt(2)/3*alpha/(alpha+1)*sin(3*phi_rad);

% Calculate the q-vector space curve in Rh coordinates and transform
% it to the experimental coordinate system
Q = R'*[sin(theta).*cos(phi_rad); sin(theta).*sin(phi_rad); cos(theta)];

% Calculate the azimuth angle in the Experimental system
chi_rad = acos(Q(3,:));
chi_deg = chi_rad * 180/pi;
% Calculate the polar angle in the Exp system
omega_rad = atan2(Q(2,:),Q(1,:));
%acos(Q(1,:)./(sqrt(Q(1,:).^2 + Q(2,:).^2)));
%
omega_deg = omega_rad*180/pi;
omega_deg_eq = linspace(0,359,360);
omega_rad_eq = omega_deg_eq*pi/180;
%

%% Plot omega dependence of the angles
Ang = [omega_deg', chi_deg'];
Ang_s = sortrows(Ang, 1);
figure(42342)
plot(Ang_s(:,1), Ang_s(:,2))
chi_deg_eq = interp1(Ang(:,1)+180, Ang(:,2), omega_deg_eq, 'linear', 'extrap');
chi_rad_eq = chi_deg_eq*pi/180;
chi_deg_eq_112 = chi_deg_eq(mod(omega_deg_eq+180,360)+1);
chi_rad_eq_112 = chi_deg_eq_112*pi/180;

figure(1)
subplot(1,4,1)
plot(omega_deg+90, phi_deg, 'k.') 
subplot(1,4,2)
plot(omega_deg+90, theta*180/pi, 'k.')
subplot(1,4,3)
plot(omega_deg+90, chi_deg, 'k.') 
hold on
plot(omega_deg_eq-180+90, chi_deg_eq, 'r-') 

subplot(1,4,4)
plot(omega_deg_eq, omega_deg, 'k.') 


figure(2)
X=sin(chi_rad_eq).*cos(omega_rad_eq+pi);
Y=sin(chi_rad_eq).*sin(omega_rad_eq+pi);
Z=cos(chi_rad_eq);
scatter3(Q(1,:), Q(2,:), Q(3,:))
hold on
scatter3(X, Y, Z)
axis square
box on

figure(5)
Q_lab = A_1m10*Q;
Q_lab2 = A_1m10*[X; Y; Z];
scatter3(Q_lab(1,:), Q_lab(2,:), Q_lab(3,:))
hold on
scatter3(Q_lab2(1,:), Q_lab2(2,:), Q_lab2(3,:))
axis square
box on

