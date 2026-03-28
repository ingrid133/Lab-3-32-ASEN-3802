% Adam Cobb
% Asen 3802
% Lab 3 Part 1 Task 1
% NACA_AIRFOILS Build boundary points 
% for a NACA 4-digit airfoil.
% 03/27/2026

clc; clear; close all;

m = 4/100;
p = 4/10; % Maximum camber position
t = 15/100; % Maximum thickness
c = 1; % Chord length

% N = panels per surface (produces N+1 chordwise nodes).
% Output length = 2*N+1.
N = 100; 


[x_b,y_b] = NACA_Airfoils(m,p,t,c,N);

% Plot the airfoil shape
figure;
plot(x_b, y_b, 'b-', 'LineWidth', 2);
axis equal;
xlabel('Chord Position (x)');
ylabel('Thickness (y)');
title('NACA Airfoil Shape');
grid on;

function [x_b, y_b] = NACA_Airfoils(m,p,t,c,N)

beta = linspace(0,pi,N+1)'; % cosine-spaced angular distribution from 0 to pi
x = (c/2) * (1-cos(beta));

xc = x/c;

yt = (t*c/0.2) .* ( 0.2969*sqrt(xc) ...
    - 0.1260*xc - 0.3516*xc.^2 + 0.2843*xc.^3 - 0.1036*xc.^4 );

% ensure closed trailing edge
yt(end) = 0;

% mean camber line yc and its derivative dyc_dx
yc = zeros(size(x));
dyc_dx = zeros(size(x));
if p > 0
    idx1 = x < p*c;
    idx2 = ~idx1;
    yc(idx1) = (m/(p^2)) .* (2*p*(x(idx1)/c) - (x(idx1)/c).^2) * c;
    yc(idx2) = (m/((1-p)^2)) .* ((1 - 2*p) + 2*p*(x(idx2)/c) - (x(idx2)/c).^2) * c;
    dyc_dx(idx1) = (2*m/(p^2)) .* (p - x(idx1)/c);
    dyc_dx(idx2) = (2*m/((1-p)^2)) .* (p - x(idx2)/c);
else
    yc(:) = 0;
    dyc_dx(:) = 0;
end

% local angle (theta)
theta = atan(dyc_dx);

% upper and lower surface coordinates
xu = x - yt .* sin(theta);
yu = yc + yt .* cos(theta);
xl = x + yt .* sin(theta);
yl = yc - yt .* cos(theta);

% assemble boundary clockwise starting at trailing edge (x = c)
% trailing edge is at the last index (x = c)
TE_x = xu(end);
TE_y = yu(end);  % equals yl(end) for perfect closed TE (may be tiny diff numerically)

% upper surface from TE -> LE: indices N+1 down to 1
upper_x = xu(end:-1:1);
upper_y = yu(end:-1:1);

% lower surface from LE -> TE: indices 2..N+1 (skip duplicate LE at index 1)
lower_x = xl(2:end);
lower_y = yl(2:end);

% combine and close at TE (repeat TE as last point)
x_b = [upper_x; lower_x; TE_x];
y_b = [upper_y; lower_y; TE_y];

% ensure column vectors
x_b = x_b(:);
y_b = y_b(:);

end
