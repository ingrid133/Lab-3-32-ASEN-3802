%% ASEN 3802 - Aerodynamics Lab - Main
% Determine effect of thickness and camber on lift 
%
% Author: Shayna Brower
% Collaborators: Adam Cobb, Luca Lungeanu, Ingrid Paska
% Date created: 2/27/26  Date revised: 4/22/26

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
%print("thickness comparison","-dpng",'-r300');

% solve for zero lift aoa with interpolation
zero_lift_aoa_vortex_0006 = interp1(cl_0006,angle_of_attack,0); % degrees
zero_lift_aoa_vortex_0012 = interp1(cl_0012_task3,angle_of_attack,0); % degrees
zero_lift_aoa_vortex_0018 = interp1(cl_0018,angle_of_attack,0); % degrees
zero_lift_aoa_experimental_0006 = interp1(data_0006(:,2),data_0006(:,1),0); % degrees
zero_lift_aoa_experimental_0012 = interp1(data_0012(:,2),data_0012(:,1),0); % degrees

% solve for lift slope per degree using polyfit
% slope for linear region (-8* <= aoa <= 8*) for experimental data
idx_exp_0006 = (data_0006(:,1) >= -8) & (data_0006(:,1) <= 8);
idx_exp_0012 = (data_0012(:,1) >= -8) & (data_0012(:,1) <= 8);

p_vortex_0006 = polyfit(angle_of_attack,cl_0006,1);
lift_slope_vortex_0006 = p_vortex_0006(1); % /degree

p_vortex_0012 = polyfit(angle_of_attack,cl_0012_task3,1); 
lift_slope_vortex_0012 = p_vortex_0012(1); % /degree

p_vortex_0018 = polyfit(angle_of_attack,cl_0018,1);
lift_slope_vortex_0018 = p_vortex_0018(1); % /degree

p_exp_0006 = polyfit(data_0006(idx_exp_0006,1),data_0006(idx_exp_0006,2),1);
lift_slope_exp_0006 = p_exp_0006(1); % /degree

p_exp_0012 = polyfit(data_0012(idx_exp_0012,1),data_0012(idx_exp_0012,2),1);
lift_slope_exp_0012 = p_exp_0012(1); % /degree

p_thin_airfoil = polyfit(angle_of_attack,cl_thin_airfoil,1);
lift_slope_thin_airfoil = p_thin_airfoil(1); % /degree


%% Task 4

NACA_2412 = [2,4,12];
NACA_4412 = [4,4,12];

% call NACA Airfoils using same N from task 2 and 3 NACA 0012 is solved in task 3
[x_b_2412,y_b_2412] = NACA_Airfoils(NACA_2412,1,N_3,'NACA 2412',0);
[x_b_4412,y_b_4412] = NACA_Airfoils(NACA_4412,1,N_3,'NACA 4412',0);

% call vortex panel for 2412 and 4412 using same range angle of attack from task 3 cl_0012 is the same as in task 3
for i=1:length(angle_of_attack)
    cl_2412(i) = Vortex_Panel(x_b_2412,y_b_2412,1,angle_of_attack(i));
    cl_4412(i) = Vortex_Panel(x_b_4412,y_b_4412,1,angle_of_attack(i));
end

% Thin airfoil theory
[zero_lift_aoa_TAT_0012,lift_slope_TAT_0012,cl_TAT_0012] = ThinAirfoilTheory(NACA_0012,1,angle_of_attack);
[zero_lift_aoa_TAT_2412,lift_slope_TAT_2412,cl_TAT_2412] = ThinAirfoilTheory(NACA_2412,1,angle_of_attack);
[zero_lift_aoa_TAT_4412,lift_slope_TAT_4412,cl_TAT_4412] = ThinAirfoilTheory(NACA_4412,1,angle_of_attack);

% Experimental data from NACA charts and digitizer
data_2412 = readmatrix("NACA2412_exp.txt");
data_4412 = readmatrix("NACA4412_exp.txt");

figure;
plot(angle_of_attack,cl_0012_task3,'LineWidth',1.4,'Color','c');
hold on;
plot(data_0012(:,1),data_0012(:,2),':o','Color','#6495ED');
plot(angle_of_attack,cl_TAT_0012,'--b','LineWidth',1.2);
plot(angle_of_attack,cl_2412,'LineWidth',1.4,'Color','m');
plot(data_2412(:,1),data_2412(:,2),':o','Color','#FF1493');
plot(angle_of_attack,cl_TAT_2412,'--r','LineWidth',1.4);
plot(angle_of_attack,cl_4412,'LineWidth',1,'Color','y');
plot(data_4412(:,1),data_4412(:,2),':o','Color',"#B8860B");
plot(angle_of_attack,cl_TAT_4412,'--','LineWidth',1.2,'Color','#FFA500');
grid on;
legend("NACA 0012 Vortex Panel","NACA 0012 Experimental",'NACA 0012 Thin Airfoil',"NACA 2412 Vortex Panel","NACA 2412 Experimental",'NACA 2412 Thin Airfoil',"NACA 4412 Vortex Panel","NACA 4412 Experimental",'NACA 4412 Thin Airfoil',"Location","southeast");
xlabel("Angle of Attack (degrees)");
ylabel("Lift Coefficient");
title("Lift Coefficient vs Angle of Attack for Airfoils of Various Amounts of Camber");
%print("camber comparison","-dpng",'-r300');

% solve for zero lift aoa with interpolation 0012 is solved in task 3
zero_lift_aoa_vortex_2412 = interp1(cl_2412,angle_of_attack,0); % degrees
zero_lift_aoa_vortex_4412 = interp1(cl_4412,angle_of_attack,0); % degrees
zero_lift_aoa_experimental_2412 = interp1(data_2412(:,2),data_2412(:,1),0); % degrees
zero_lift_aoa_experimental_4412 = interp1(data_4412(1:16,2),data_4412(1:16,1),0); % degrees

% solve for lift slope per degree using polyfit
% slope for linear region (-8* <= aoa <= 8*) for experimental data
idx_exp_2412 = (data_2412(:,1) >= -8) & (data_2412(:,1) <= 8);
idx_exp_4412 = (data_4412(:,1) >= -8) & (data_4412(:,1) <= 8);

p_vortex_2412 = polyfit(angle_of_attack,cl_2412,1);
lift_slope_vortex_2412 = p_vortex_2412(1); % /degree

p_vortex_4412 = polyfit(angle_of_attack,cl_4412,1);
lift_slope_vortex_4412 = p_vortex_4412(1); % /degree

p_exp_2412 = polyfit(data_2412(idx_exp_2412,1),data_2412(idx_exp_2412,2),1);
lift_slope_exp_2412 = p_exp_2412(1); % /degree

p_exp_4412 = polyfit(data_4412(idx_exp_4412,1),data_4412(idx_exp_4412,2),1);
lift_slope_exp_4412 = p_exp_4412(1); % /degree


%% Part 2 Task 1

N_part2 = 50;
b = 10; %ft
a0_r = 2*pi; % /rad
a0_t = 2*pi; % /rad
aero_r = 0;  % deg
aero_t = 0;  % deg
geo_r = 5;   % deg
geo_t = 5;   % deg

AR = [4,6,8,10]; % AR values from fig 5.20
taper_ratio = linspace(0,1,1000);


% use for loops to solve for ct/cr vs induced drag factor for each AR
for i=1:length(AR)
    Aspect_ratio = AR(i);

    for j=1:length(taper_ratio)
        ct_cr = taper_ratio(j);
        
        % solve for c_r and c_t to input into fuction 
        c_r = 2*b/(Aspect_ratio*(1+ct_cr));
        c_t = ct_cr*c_r;
        
        % call PLLT function for each c_r and c_t
        [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N_part2);
        
        % solve for induced drag factor
        delta(j,i) = (1/e) - 1;
    end
end

% reproduce figure 5.20
figure;
for i=1:length(AR)
    plot(taper_ratio,delta(:,i),'LineWidth',1.5);
    hold on;
end
xlabel('Taper ratio, ct/cr');
ylabel('\delta');
legend('AR = 4', 'AR = 6','AR = 8', 'AR = 10',"Location",'best');
%print('5_20 plot','-dpng','-r300');

%% Part 3

% givens
aoa_3 = 4; % degrees
b_3 = 400/12; %ft
cr_3 = 64/12; %ft
ct_3 = 44.5/12; %ft
geo_r_3 = aoa_3 + 1; %degree 
geo_t_3 = aoa_3 + 0; %degree 
a0_r_3 = (lift_slope_vortex_2412*180)/pi; % /rad NACA 2412
a0_t_3 = (lift_slope_vortex_0012*180)/pi; % /rad NACA 0012
aero_r_3 = zero_lift_aoa_vortex_2412;  % deg NACA 2412
aero_t_3 = zero_lift_aoa_vortex_0012;  % deg NACA 0012

N_ref_part_3 = 600;
N_part3 = 1:1:N_ref_part_3;

% solve for CL and CDi at a large N for reference for convergence study
[e_3_ref_3,c_L_ref_3,c_Di_ref_3] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,N_ref_part_3);

% solve for CL error values
CL_10_percent = 0.1*c_L_ref_3;
CL_1_percent = 0.01*c_L_ref_3;
CL_1_tenth_percent = 0.001*c_L_ref_3;

% put error values in a vector
CL_error_values = [CL_10_percent,CL_1_percent,CL_1_tenth_percent];

% find number of terms for 10% error, 1% error and .1% error in CL
for j = 1:length(CL_error_values)
    a = 0; %restart counter before for loop
    for i= 1:length(N_part3)
    a = a+1;
    [e_3_part3(a),c_L_part3(a),c_Di_part3(a)] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,N_part3(i));
    N_values_part_3_CL(a) = i; % store N values for plotting
    if (c_L_part3(a) >= (c_L_ref_3 - CL_error_values(j))) && (c_L_part3(a) <= (c_L_ref_3 + CL_error_values(j)))
        min_odd_terms_CL(j) = N_part3(i); % N for 10%, 1%, and .1% error
        CL_at_min_terms(j) = c_L_part3(a); % CL for 10%, 1%, and .1% error 
        break;
    end
    end
end

% solve for the remaining CL terms until reference N_ref is reached iterating N by 10
for i = min_odd_terms_CL(3) + 10:10:N_ref_part_3
    a = a+1;
    [e_3_part3(a),c_L_part3(a),c_Di_part3(a)] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,N_part3(i));
    N_values_part_3_CL(a) = i;
end


% solve for CDi error values
CDi_10_percent = 0.1*c_Di_ref_3;
CDi_1_percent = 0.01*c_Di_ref_3;
CDi_1_tenth_percent = 0.001*c_Di_ref_3;

CDi_error_values = [CDi_10_percent,CDi_1_percent,CDi_1_tenth_percent];

% find number of terms for 10%, 1%, and .1% error in CDi
for j = 1:length(CDi_error_values)
    a = 0; %restart counter before for loop
    for i= 1:length(N_part3)
        a = a+1;
        [e_3_part3(a),c_L_part3_CDi(a),c_Di_part3(a)] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,N_part3(i));
        N_values_part_3_CDi(a) = i; %store N values for plotting
        if (c_Di_part3(a) >= (c_Di_ref_3 - CDi_error_values(j))) && (c_Di_part3(a) <= (c_Di_ref_3 + CDi_error_values(j)))
            min_odd_terms_CDi(j) = N_part3(i); % N for 10%, 1%, and .1% error
            CDi_at_min_terms(j) = c_Di_part3(a); % CDi for 10%, 1%, and .1% error
            break;
        end
    end
end


% solve for the remaining CDi terms until reference N_ref is reached iterating N by 10
for i = min_odd_terms_CDi(3) + 10:10:N_ref_part_3
    a = a+1;
    [e_3_part3(a),c_L_part3_CDi(a),c_Di_part3(a)] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,N_part3(i));
    N_values_part_3_CDi(a) = i;
end

% plots for deliverable 2
figure;
plot(N_values_part_3_CL,c_L_part3,'LineWidth',1.5);
hold on;
xline(min_odd_terms_CL(1),'r');
xline(min_odd_terms_CL(2),'Color','#33CC33');
xline(min_odd_terms_CL(3),'k');
xlabel('Number of Odd Terms');
ylabel('Coefficient of Lift');
title("CL vs. Number of Odd Terms");
legend('CL','10% error','1% error','0.1% error');
xlim([0 200]);
%print("CL","-dpng",'-r300');

figure;
plot(N_values_part_3_CDi,c_Di_part3,'LineWidth',1.5);
hold on;
xline(min_odd_terms_CDi(1),'r');
xline(min_odd_terms_CDi(2),'Color','#33CC33');
xline(min_odd_terms_CDi(3),'k');
xlabel('Number of Odd Terms');
ylabel('Induced Drag Coefficient');
title("CDi vs. Number of Odd Terms");
legend('CDi','10% error','1% error','0.1% error');
xlim([0 200]);
%print("CDi","-dpng",'-r300');

%% Deliverable 3

% call PLLT with N for .1% error CDi
[e_del3,c_L_del3,c_Di_del3] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_3,geo_r_3,min_odd_terms_CDi(3));

V = 100*1.68781; %ft/s
h = 10000; % ft
rho = 1.7556e-3; % slugs/ft^3
S = b_3*(cr_3 + ct_3)/2; % ft^2
q = 0.5*rho*(V^2); % lb/ft^2

% solve for Lift and induced drag
L = q*S*c_L_del3; % lbf
Di = q*S*c_Di_del3; % lbf

% solve for cd at aoa = 4 degrees using NACA chart data
cl_data_0012 = readmatrix('NACA0012_exp.txt');
cd_data_0012 = readmatrix('NACA0012_cd_data.txt');
cl_data_2412 = readmatrix('NACA2412_exp.txt');
cd_data_2412 = readmatrix('NACA2412_cd_data.txt');
cl_4_0012 = interp1(cl_data_0012(:,1),cl_data_0012(:,2),4);
cl_4_2412 = interp1(cl_data_2412(:,1),cl_data_2412(:,2),4);
cd_4_0012 = interp1(cd_data_0012(:,1),cd_data_0012(:,2),cl_4_0012);
cd_4_2412 = interp1(cd_data_2412(:,1),cd_data_2412(:,2),cl_4_2412);

% solve for average section drag coefficient given 2412 at the root and 0012 at the tip to solve for total drag 
cd_del3_avg = (cr_3*cd_4_2412 + ct_3*cd_4_0012) / (cr_3 + ct_3);

% solve for total L/D
CD = c_Di_del3 + cd_del3_avg;
D = CD*S*q; %lbf
L_D_efficiency = L/D;

%% Deliverable 4

aoa_del4 = -8:1:8;

% use for loop to solve for CD, CDi, and Cd as angle of attack ranges from -8 to 8 degrees
for i=1:length(aoa_del4)
    cl_del4_0012(i) = interp1(cl_data_0012(:,1),cl_data_0012(:,2),aoa_del4(i));
    cl_del4_2412(i) = interp1(cl_data_2412(:,1),cl_data_2412(:,2),aoa_del4(i));
    cd_del4_0012(i) = interp1(cd_data_0012(:,1),cd_data_0012(:,2),cl_del4_0012(i));
    cd_del4_2412(i) = interp1(cd_data_2412(:,1),cd_data_2412(:,2),cl_del4_2412(i));
    
    % update geometric angle of attack as aoa changes
    geo_r_del4(i) = aoa_del4(i) + 1; %degree 
    geo_t_del4(i) = aoa_del4(i) + 0; %degree 
    
    % solve for CDi using PLLT, average cd, and total CD as angle of attack changes
    [e_del4(i),c_L_del4(i),c_Di_del4(i)] = PLLT(b_3,a0_t_3,a0_r_3,ct_3,cr_3,aero_t_3,aero_r_3,geo_t_del4(i),geo_r_del4(i),min_odd_terms_CDi(3));
    cd_del4_avg(i) = (cr_3*cd_del4_2412(i) + ct_3*cd_del4_0012(i)) / (cr_3 + ct_3);
    CD_del4(i) = c_Di_del4(i) + cd_del4_avg(i);
end

% plot for deliverable 4 
figure;
plot(aoa_del4,cd_del4_avg,'LineWidth',1.5);
hold on;
plot(aoa_del4,c_Di_del4,'LineWidth',1.5);
plot(aoa_del4,CD_del4,'LineWidth',1.5);
xlabel('Angle of Attack (degrees)');
ylabel('Drag Coefficient');
legend("Profile drag coefficient",'Induced drag coefficent','Total drag coefficient');
title('Drag Coefficient vs Angle of Attack');
print("CD vs aoa","-dpng",'-r300');


%% Deliverable 5

% solve for L,D and L/D
L_del5 = c_L_del4.*q.*S;
D_del5 = CD_del4.*q.*S;
L_D_efficiency_del5 = L_del5./D_del5;

% plot for deliverable 5
figure; 
plot(aoa_del4,L_D_efficiency_del5,"LineWidth",1.5);
hold on;
xlabel('Angle of Attack (degrees)');
ylabel('Aerodynamic Efficency (L/D)');
title('L/D vs Angle of Attack');
print("L_D vs aoa","-dpng","-r300");