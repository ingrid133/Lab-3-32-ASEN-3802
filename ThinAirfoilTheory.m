function [zero_lift_angle_of_attack,lift_slope_per_degree,cl] = ThinAirfoilTheory(NACA, c, alpha_deg)
%THINAIRFOILTHEORY solve for cl given angle of attack using thin airfoil theory
%
% Author: Shayna Brower
% Collaborators: Adam Cobb, Luca Lungeanu, Ingrid Paska
% Date created: 4/8/26  Date revised: 4/8/26

m = NACA(1)/100;
p = NACA(2)/10;

alpha_rad = deg2rad(alpha_deg);

% symmetric airfoil
if (m == 0) && (p ==0)
    zero_lift_angle_of_attack = 0; %degrees
    cl = 2*pi*alpha_rad;
    P = polyfit(alpha_deg,cl,1);
    lift_slope_per_degree = P(1);
    return;
end


% cambered airfoil
theta = linspace(0,pi,2000);
x = (c/2) * (1-cos(theta));

for i = 1:length(x)
    if (p ~= 0)
        if (x(i) >= 0) && (x(i) < (p*c))
            dyc_dx(i) = ((2*m)/(p^2))*(p - (x(i)/c));
        else
            dyc_dx(i) = ((2*m)/(1-p)^2)*(p - (x(i)/c));
        end
    end
end

% Thin airfoil theory zero-lift angle in radians
zero_lift_angle_of_attack_rad = -(1/pi) * trapz(theta, dyc_dx .* (cos(theta) - 1));

% Convert to degrees
zero_lift_angle_of_attack = rad2deg(zero_lift_angle_of_attack_rad); %degrees

% Build lift curve
cl = 2*pi*(alpha_rad - zero_lift_angle_of_attack_rad);

% lift slope
P = polyfit(alpha_deg,cl,1);
lift_slope_per_degree = P(1);

end