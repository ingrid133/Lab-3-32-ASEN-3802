function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%e is the span efficiency factor (to be computed and returned)
%C_L is the coefficient of lift (to be computed and returned)
%C_Di is the induced coefficient of drag (to be computed and returned)
%b is the span (in feet
%a0_t is the cross-sectional lift slope at the tips (per radian)
%a0_r is the cross-sectional lift slope at the root (per radian)
%c_t is the chord at the tips (in feet)
%c_r is the chord at the root (in feet)
%aero_t is the zero-lift angle of attack at the tips (in degrees)
%aero_r is the zero-lift angle of attack at the root (in degrees)
%geo_t is the geometric angle of attack at the tips (in degrees)
%geo_r is the geometric angle of attack at the root (in degrees)
%N is the number of odd terms to include in the series expansion for circulation

%calc terms
S= b*(c_r + c_t)/2;
AR= b^2/S;
A_matrix = zeros(N,N);
RHS = zeros(N,1);
%need be in rad or will interp wrong
aero_t= (pi/180)*aero_t;
aero_r= (pi/180)*aero_r;
geo_t = (pi/180)*geo_t;
geo_r=(pi/180)*geo_r;

%theta loop for the diff theta
theta = zeros(N,1);
for i= 1:N
    theta(i)= (i*pi)/(2*N);
end

for i= 1:N
    frac = cos(theta(i));
    c=c_r +(c_t-c_r)*frac;
    a0=a0_r +(a0_t- a0_r)*frac;
    alp_L0= aero_r+ (aero_t-aero_r)*frac;
    alp = geo_r +(geo_t-geo_r)*frac;
    RHS(i)= alp-alp_L0;
    for j=1:N
        n=2*j -1;
        A_matrix(i,j)= (4*b/(a0*c))*sin(n*theta(i))+ n*sin(theta(i)*n)/sin(theta(i));

    end

end
A = A_matrix\RHS;
c_L= pi*AR*A(1);
sum=0;
for j =1:N
    n=2*j -1;
    sum = sum +n*(A(j)^2);
end
c_Di =pi*AR*sum;
e =(c_L^2)/(pi*AR*c_Di);
end