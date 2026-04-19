function [x_b,y_b] = NACA_Airfoils(NACA,c,N,title_name,plot_on_off)
% NACA_AIRFOILS create a plot of boundary points for NACA airfoils
% Using vortex panel method solve for and plot the x and y locations of 
% the boundary points for a given 4-digit NACA airfoil
%
% Author: Shayna Brower
% Collaborators: Adam Cobb, Luca Lungeanu, Ingrid Paska
% Date created: 2/27/26  Date revised: 4/6/26

% Task 1

% solve for m, p, and t given an NACA airfoil 
m = NACA(1)/100;
p = NACA(2)/10;
t = NACA(3)/100;

% use vortex panel method for x locations
beta = linspace(0,pi,N);
x = (c/2)*(1 - cos(beta));

% solve for thickness distribution 
y_t = ((t*c)/.2) * ( (.2969.*sqrt(x./c)) - (.1260.*(x./c)) - (.3516.*((x./c).^2)) + (.2843.*((x./c).^3))  - (.1036.*((x./c).^4)));

% solve for y values for mean camber line and derivative of camber line with respect to x if given a cambered airfoil
y_c = zeros(1,N);
dyc_dx = zeros(1,N);
for i=1:N
    if (p ~= 0)
        if (x(i) >= 0) && (x(i) < (p*c))
            y_c(i) = m*(x(i)/(p^2))*((2*p) - (x(i)/c));
            dyc_dx(i) = ((2*m)/(p^2))*(p - (x(i)/c));
        else
            y_c(i) = m*((c-x(i))/((1-p)^2))*( 1 + (x(i)/c) - (2*p));
            dyc_dx(i) = ((2*m)/(1-p)^2)*(p - (x(i)/c));
        end
    end
end

% solve for local angle
local_angle = atan(dyc_dx);

% solve for x and y upper and lower surfaces of the airfoil
x_upper = x - y_t.*sin(local_angle);
x_lower = x + y_t.*sin(local_angle);

y_upper = y_c + y_t.*cos(local_angle);
y_lower = y_c - y_t.*cos(local_angle);


% create vectors for the x and y values of the boundary points clockwise starting from the trailing edge
x_b = zeros(1,2*N-1);
y_b = zeros(1,2*N-1);
for i = 1:N
    x_b(i) = x_lower((N+1)-i); 
    y_b(i) = y_lower((N+1)-i); 
end
for i = 2:N
    x_b(N+i-1) = x_upper(i);
    y_b(N+i-1) = y_upper(i);
end

% plot airfoil
if plot_on_off == 1
    figure;
    plot(x_b,y_b,'-*','LineWidth',1.5);
    hold on;
    axis equal;
    grid on;
    if (m == 0)
        legend('Airfoil Boundary Points', 'Location','best');
    else
        plot(x,y_c,'LineWidth',1.5);
        legend('Airfoil Boundary Points','Mean Camber Line', 'Location','best');
    end
    title(title_name);
    xlabel('Chord Position');
    ylabel('Thickness');
end

end