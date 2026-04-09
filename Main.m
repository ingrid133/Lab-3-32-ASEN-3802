%% ASEN 3802 - Aerodynamics Lab - Main
% 
% Author: Shayna Brower
% Collaborators: Adam Cobb, Luca Lungeanu, Ingrid Paska
% Date created: 2/27/26  Date revised: 4/6/26

clear;
clc;
close all;

%% Task 1

NACA_0021 = [0,0,21];
NACA_2421 = [2,4,21];
N = 51; % 50 panels for each surface

% call function for each airfoil
[x_b_0021,y_b_0021] = NACA_Airfoils(NACA_0021,1,N,'NACA 0021',1);
%print("NACA_0021","-dpng","-r300");

[x_b_2421,y_b_2421] = NACA_Airfoils(NACA_2421,1,N,'NACA 2421',1);
%print("NACA_2421","-dpng","-r300");

%% Task 2

N_ref = 601;

% givens
NACA_0012 = [0,0,12];
alpha = 12; %degrees

% call function for "exact" solution using very large number of panels
[x_b_0012_ref,y_b_0012_ref] = NACA_Airfoils(NACA_0012,1,N_ref,'NACA 0012',0);

% call vortex panel function
cl_0012_ref = Vortex_Panel(x_b_0012_ref,y_b_0012_ref,1,alpha);
disp(['Exact cl for NACA airfoil 0012 = ', num2str(cl_0012_ref)]);

% solve for 1% relative error in cl 
cl_max_0012 = cl_0012_ref + (.01*cl_0012_ref);
cl_min_0012 = cl_0012_ref - (.01*cl_0012_ref);

% solve for the minimum number of panels that is within the +- 1% relative error range
a = 0;
for i = 5:N_ref
    a = a+1;
    [x_b_0012_temp,y_b_0012_temp] = NACA_Airfoils(NACA_0012,1,i,'NACA 0012',0);
    cl_0012(a) = Vortex_Panel(x_b_0012_temp,y_b_0012_temp,1,alpha);
    N_values(a) = i;
    if (cl_0012(a) >= cl_min_0012) && (cl_0012(a) <= cl_max_0012)
        min_number_of_boundary_points_per_surface = i;
        cl_approx = cl_0012(a);
        disp(['Minimum total number of panels = ', num2str(2*(min_number_of_boundary_points_per_surface-1))]);
        break;
    end
end

% solve for cl_values after min panels iterating boundary points by 20 every time
for i = min_number_of_boundary_points_per_surface + 20:20:N_ref
    a = a+1;
    [x_b_0012_temp,y_b_0012_temp] = NACA_Airfoils(NACA_0012,1,i,'NACA 0012',0);
    cl_0012(a) = Vortex_Panel(x_b_0012_temp,y_b_0012_temp,1,alpha);
    N_values(a) = i;
end

% solve for total number of panels as N increases
number_of_panels = 2*(N_values-1);

% convergence plot
figure;
plot(number_of_panels,cl_0012,'LineWidth',1.5,'Color','c');
hold on;
plot((2*(N_ref-1)),cl_0012_ref,'*m');
xline(2*(min_number_of_boundary_points_per_surface-1),'--y','LineWidth',1.5);
yline(cl_min_0012,'--r','LineWidth',1.5);
yline(cl_max_0012,'--r','LineWidth',1.5);
grid on;
xlim([0 1250]);
legend('Lift coefficient vs total number of panels','Exact solution using large number of panels','Minimum number of panels','1% error bounds','Location','southeast');
xlabel("Total Number of Panels");
ylabel("Lift Coefficient");
title("Convergence of Predicted Lift Coefficient with Respect to Total Number of Panels");
%print("convergence plot","-dpng","-r300");


%% Task 3

NACA_0006 = [0,0,6];
NACA_0012 = [0,0,12];
NACA_0018 = [0,0,18];

% use number of boundary points found in task 2
N_3 = min_number_of_boundary_points_per_surface; 

% create a vector for angle of attack values
angle_of_attack = -16:1:20; %degrees

[x_b_0006,y_b_0006] = NACA_Airfoils(NACA_0006,1,N_3,'NACA 0006',0);
[x_b_0012,y_b_0012] = NACA_Airfoils(NACA_0012,1,N_3,'NACA 0012',0);
[x_b_0018,y_b_0018] = NACA_Airfoils(NACA_0018,1,N_3,'NACA 0018',0);

% solve for cl with respect to alpha for all three airfoils
for i=1:length(angle_of_attack)
    cl_0006(i) = Vortex_Panel(x_b_0006,y_b_0006,1,angle_of_attack(i));
    cl_0012_task3(i) = Vortex_Panel(x_b_0012,y_b_0012,1,angle_of_attack(i));
    cl_0018(i) = Vortex_Panel(x_b_0018,y_b_0018,1,angle_of_attack(i));
end


% Thin airfoil theory
% dz/dx = 0 since all three are not cambered airfoils, so cl = 2*pi*angle of attack
% cl slope is the same for all three airfoils 
angle_of_attack_radians = (pi.*angle_of_attack)./180; % convert to radians
cl_thin_airfoil = 2*pi.*angle_of_attack_radians;

% experimental data from NACA charts and digitizer
data_0006 = readmatrix("NACA0006_exp.txt");
data_0012 = readmatrix("NACA0012_exp.txt");


% plot cl vs aoa for all three methods
figure;
plot(angle_of_attack,cl_0006,'LineWidth',1.2,'Color','m');
hold on;
plot(data_0006(:,1),data_0006(:,2),'--om');
plot(angle_of_attack,cl_0012_task3,'LineWidth',1.2,'Color','c');
plot(data_0012(:,1),data_0012(:,2),'--oc');
plot(angle_of_attack,cl_0018,'LineWidth',1,'Color','y');
plot(angle_of_attack,cl_thin_airfoil,'--w','LineWidth',1.2);
grid on;
legend("NACA 0006 Vortex Panel","NACA 0006 Experimental","NACA 0012 Vortex Panel","NACA 0012 Experimental","NACA 0018 Vortex Panel","Thin Airfoil Theory","Location","southeast");
xlabel("Angle of Attack (degrees)");
ylabel("Lift Coefficient");
title("Lift Coefficient vs Angle of Attack for Airfoils of Various Thickness");
print("thickness comparison","-dpng",'-r300');



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
    plot(x_b,y_b,'LineWidth',1.5);
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