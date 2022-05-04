
function [theta3,theta4] = fourbarpos(len,theta1,theta2,delta)
delta = -1;
    K1 = len(4)/len(1);
    K2 = len(4)/len(3);
    K3 = (len(1)^2 - len(2)^2 + len(3)^2 + len(4)^2)/(2*len(1)*len(3));
    K4 = len(4)/len(2);
    K5 = (len(3)^2 - len(4)^2 - len(1)^2 - len(2)^2)/(2*len(1)*len(2));

    A = cosd(theta2) - K1 - K2*cosd(theta2) + K3;
    B = -2*sind(theta2);
    C = K1 - (K2+1)*cosd(theta2) + K3;
    D = cosd(theta2) - K1 + K4*cosd(theta2) + K5;
    E = -2*sind(theta2);
    F = K1 + (K4-1)*cosd(theta2) + K5;

    theta4 = 2*atan2d((2*A),-B  + delta*sqrt(B^2-4*A*C));
    theta3 = 2*atan2d((2*D),-E  + delta*sqrt(E^2-4*D*F));
end