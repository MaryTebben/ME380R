% clc; clear; close all;

L1x = 79.09;
L1y = 197.68;
L1 = sqrt(197.68^2+79.09^2);
L2 = 108.86;
L3 = 89.92;
L4 = 15;
L5 = 121.64;

th2L1 = abs(atand(197.68/79.07));
theta1 = 180-th2L1;
theta2 = 90;
% one rotation every 5 seconds (set by how much of the song is chosen to be
% played)
omegar = 2*pi/5; 
alphar = 0;

%% Find cam data

% picture = ('cam picture.png');
% rgbImage = imread(picture);
% [rows, columns, numberOfColorBands] = size(rgbImage);
% 
% subplot(2, 2, 1);
% imshow(rgbImage);
% axis on;
% title('Original Color Image');
% 
% % Extract the individual red, green, and blue color channels.
% redChannel = rgbImage(:, :, 1);
% % greenChannel = rgbImage(:, :, 2);
% % blueChannel = rgbImage(:, :, 3);
% 
% % Get the binaryImage
% binaryImage = redChannel < 200;
% % Display the binary image.
% subplot(2, 2, 2);
% imshow(binaryImage);
% axis on;
% title('Binary Image');
% 
% [y,x]=find(binaryImage==1); 
% xy=[x,y]; 
% media=mean(xy); 
% subplot(2, 2, 3);
% %  imshow(binaryImage); 
% plot(x,y)
% axis ij;
% hold on 
% plot(media(:,1),media(:,2),'*b');
% hold off
% 
% % grabit is a function found on MathWorks to allow manual selection of
% % desired image coordinates
% % grabit('cam picture.png')
% subplot(2, 2, 4)
load('cam_points.mat');
load('center.mat');
cam_center = mean(center);
% plot(cam_points(:,1),cam_points(:,2))
% hold on
% plot(cam_center(:,1),cam_center(:,2),'*b')
% hold off
% axis ij;

%% Kinematic analysis

% r is cam radius starting from desired location
r = sqrt((cam_points(:,1) - cam_center(:,1)).^2 + (cam_points(:,2) - cam_center(:,2)).^2);
% from solidworks
r_start = 34.53;
scale = r_start/r(1);
r = r*scale;
crank_angle = rad2deg(linspace(0,2*pi,length(r)));
% 
% for i = 1:length(r)
%     if i == 1
%         r_dm1 = r(length(r));
%         r_dp1 = r(i+1);
%         ca_dm1 = crank_angle(length(r));
%         ca_dp1 = crank_angle(i+1);
%     elseif i == length(r)
%         r_dm1 = r(i-1);
%         r_dp1 = r(1);
%         ca_dm1 = crank_angle(i-1);
%         ca_cp1 = crank_angle(1);
%     else
%         r_dm1 = r(i-1);
%         r_dp1 = r(i+1);
%         ca_dm1 = crank_angle(i-1);
%         ca_dp1 = crank_angle(i+1);
%     end
%     r_d(i) = (r_dp1 - r_dm1)/(ca_dp1 - ca_dm1);
% end
% for i = 1:length(r_d)
%     if i == 1
%         rd_dm1 = r_d(length(r_d));
%         rd_dp1 = r_d(i+1);
%         ca_dm1 = crank_angle(length(r));
%         ca_dp1 = crank_angle(i+1);
%     elseif i == length(r_d)
%         rd_dm1 = r_d(i-1);
%         rd_dp1 = r_d(1);
%         ca_dm1 = crank_angle(i-1);
%         ca_cp1 = crank_angle(1);
%     else
%         rd_dm1 = r_d(i-1);
%         rd_dp1 = r_d(i+1);
%         ca_dm1 = crank_angle(i-1);
%         ca_dp1 = crank_angle(i+1);
%     end
%     r_dd(i) = (rd_dp1 - rd_dm1)/(ca_dp1 - ca_dm1);
% end
% 
% for i = 1:length(r)
%     syms t3 t4
%     position1 = 0 == (r(i)+L2)*cosd(theta2) - L3*cosd(t3) - L4*cosd(t4) - L1*cosd(theta1);
%     position2 = 0 == (r(i)+L2)*sind(theta2) - L3*sind(t3) - L4*sind(t4) - L1*sind(theta1);
%     
%     pos_solution = vpasolve([position1,position2],[t3,t4]);
%     t3 = double(pos_solution.t3);
%     t4 = double(pos_solution.t4);
%     theta3(i) = t3;
%     theta4(i) = t4;
% 
%     syms o3 o4
%     velocity1 = 0 == r_d(i)*cosd(theta2) + L3*o3*sind(theta3(i)) + L4*o4*sind(theta4(i));
%     velocity2 = 0 == r_d(i)*sind(theta2) - L3*o3*cosd(theta3(i)) - L4*o4*cosd(theta4(i));
% 
%     vel_solution = vpasolve([velocity1,velocity2],[o3,o4]);
%     o3 = double(vel_solution.o3);
%     o4 = double(vel_solution.o4);
%     omega3(i) = o3;
%     omega4(i) = o4;
% 
%     syms a3 a4
%     acceleration1 = 0 == r_dd(i)*cosd(theta2) + L3*omega3(i)^2*cosd(theta3(i)) + L3*a3*sind(theta3(i)) + L4*omega4(i)^2*cosd(theta4(i)) + L4*a4*sind(theta4(i));
%     acceleration2 = 0 == r_dd(i)*sind(theta2) + L3*omega3(i)^2*sind(theta3(i)) - L3*a3*cosd(theta3(i)) + L4*omega4(i)^2*sind(theta4(i)) - L4*a4*cosd(theta4(i));
% 
%     acc_solution = vpasolve([acceleration1,acceleration2],[a3,a4]);
%     a3 = double(acc_solution.a3);
%     a4 = double(acc_solution.a4);
%     alpha3(i) = a3;
%     alpha4(i) = a4;
% end

omega5 = -omega4;
alpha5 = -alpha4;

crank_angle = rad2deg(linspace(0,2*pi,length(r)));

figure(99)
plot(crank_angle,theta3,crank_angle,theta4)
legend('\theta_3','\theta_4')

P5x = L1*cosd(theta1) + L5.*cosd(-theta4);
P5y = L1*sind(theta1) + L5.*sind(-theta4);
P5 = abs(P5x + 1i*P5y);
ang5 = rad2deg(angle(P5x+1i*P5y));
V5 = L5.*omega5;
A5_n = L5.*omega5.^2;
A5_t = L5*alpha5;
A5 = A5_t + A5_n;

figure(1)
plot(P5x,P5y)%,crank_angle,ang5);
title('Point 5 Position');
axis ij; set(gca,'XDir','reverse');
xlabel('xpos'); ylabel('ypos');
figure(2)
plot(crank_angle,V5);
title('Point 5 Velocity')
xlabel('Cam rotation angle'); ylabel('Velocity (mm/s)');
figure(3)
plot(crank_angle,A5)
title('Point 5 Acceleration')
xlabel('Cam rotation angle'); ylabel('Acceleration (mm/s^2)');