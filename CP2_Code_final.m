%CP2
clc;
clear all;
close all;

%  ------------------------------
%  Constants
%  ------------------------------
kp = 1.0;     % kappa
h  = 1.0;
Ti = 5.0;
Q1 = 10.0;
Q2_ref = 5.0;
n = 12;

% ----------------------------------------
% Verification for Thomas algorithm
% ----------------------------------------
A =[2 -1  0  0 ;1  2 -1  0. ;0  1  2 -1 ;0  0  1  2];
d = [0.0 6.0 2.0 -5.0];
[a, b, c] = LUdecomp(A);
[x] = LUsolve(a, b, c, d);


% ----------------------------------------
% (4)(1)  This is how the heatSolver is used.
% ----------------------------------------
dts = [0.01 0.05 0.25];
[y1, t1, u1] = heatSolver(0.1, dts(1), Q2_ref); %using the implimented Heatsolver
[y2, t2, u2] = heatSolver(0.1, dts(2), Q2_ref);%using the implimented Heatsolver
[y3, t3, u3] = heatSolver(0.1, dts(3), Q2_ref);%using the implimented Heatsolver


% ----------------------------------------
% (4)(2) For you to figure out
% ----------------------------------------
Ns=[5,10,20,30,60,120,200];
for i=1:length(Ns)
    t=Ns(i)*0.01;
    plot(y1,anaSol(y1,t,n,Q2_ref),"--",'DisplayName',['Analyticalsol for Time ' num2str(t)])
    hold on
    plot(y1,u1(Ns(i),:) ,'DisplayName',['Numericsol for Time ' num2str(t)])

end
xlabel("y");
ylabel("temp");
legend();
%--
%the plot above shows very little diviation from the analytical to the
%numerical solution and it can be found that as the time increases the
% temperature also increases throughout the beam
    
% ----------------------------------------
% (4)(3) For you to figure out
% ----------------------------------------
figure()
subplot(2,1,1);
plot(t1,u1(:,1),'-','DisplayName',['Timestep ' num2str(0.01)]);
hold on
plot(t2,u2(:,1),'-','DisplayName',['Timestep  ' num2str(0.05)]);
plot(t3,u3(:,1),'DisplayName',['Timestep ' num2str(0.25)]);
plot(t1,anaSol(0,t1,n,Q2_ref),'--','DisplayName',['Timestep_analytical ' num2str(0.01)]);
plot(t2,anaSol(0,t2,n,Q2_ref),'--','DisplayName',['Timestep_analytical ' num2str(0.05)]);
plot(t3,anaSol(0,t3,n,Q2_ref),'--','DisplayName',['Timestep_analytical ' num2str(0.25)]);
ylabel("temp")
legend();

subplot(2,1,2)
plot(t1,u1(:,1)-anaSol(0,t1,n,Q2_ref),'-');
hold on
plot(t2,u2(:,1)-anaSol(0,t1,n,Q2_ref),'-');
plot(t3,u3(:,1)-anaSol(0,t1,n,Q2_ref),'-');
%from the above graph it can be observed that as the delta time increases
%the more the deviasion increases from the numerical and analytical curve
%indicating that the smaller time step is much more accurate in this
%analysis.

% ----------------------------------------
% (4)(4) For you to figure out
% ----------------------------------------

Qs = linspace(5,15,31);
us=[];
for i=1:length(Qs)
    [~,~,u]= heatSolver(0.1,0.5,Qs(i));
    us=[us;u(length(u),1)];
end

figure()
plot(Qs,us,'b-','DisplayName',"Numerical Solution");
xlabel("Q");
ylabel("Max Temp");
% from the above graph we can see that the  when umax=5 the corresponding Q
% value is 10 and based on these conditions with q2 =10, we can design a
% new system to achive these conditions.

function [F, T, S, B1, B2] = initProblem(Ny, dy, dt, Q2)
%     """
%     A helper function called by *heatSolver* that prepares the quantities in Eq. (11)
%     Input:
%     Ny: Number of steps in y-direction, an integer
%     dy: Step size in y-direction, a float
%     dt: Time step size, a float
%     Q2: Heat flux of the cooling system, default Q2=5.0
%     Output:
%     F: The forcing vector, (_Ny-1)-element vector
%     T: The LHS matrix, (_Ny-1)x(_Ny-1) matrix
%     S: The RHS matrix, (_Ny-1)x(_Ny-1) matrix
%     B1: Coefficients for computing u(0,t), 3-element vector
%     B2: Coefficients for computing u(1,t), 3-element vector
%     """

Q1=10;
nu=1;
N=Ny-1;
f=zeros(N,1);
f(1)=Q1*2/3;
f(N)=-Q2*2/3;
A= -2*diag(ones(1,N),0)+diag(ones(1,N-1),1)+diag(ones(1,N-1),-1);
A0=zeros(1,N);
A0(1)=-2/3;
A0(2)=-A0(1);
A(1,:)=A0;
A(N,:)=fliplr(A0);
I=diag(ones(1,N),0);
r=(nu/2)*(dt/dy^2);
F=dt/dy*nu*f;
T=I-r*A;
S=I+r*A;

% ----------------------------------------
% TODO: Construct the F, T, S arrays
% ----------------------------------------
% F = ???
% T = ???
% S = ???

% B vectors - already provided here
B1 =[Q1*(2*dy)/3.0 4.0/3.0 -1.0/3.0];
B2 =[-Q2*(2*dy)/3.0 4.0/3.0 -1.0/3.0];
end

function [a, b, c] = LUdecomp(T)
%     LU decomposition of a tridiagonal matrix in the form of three arrays.
%     Input:
%     T:  Tridiagonal matrix to decompose, NxN matrix
%     Output:
%     a: Main diagonal of U matrix, N-element array
%     b: Lower diagonal of L matrix, (N-1)-element array
%     c: Upper diagonal of U matrix, (N-1)-element array


a = diag(T);      % Initialize a by main diagonal of T
b = diag(T, -1);  % Initialize b by lower diagonal of T
c = diag(T, 1);   % Initialize c by upper diagonal of T
N = length(a);
for i=2:N
    % ----------------------------------------
    % TODO: Complete the loop to compute arrays of a and b
    % ----------------------------------------
    b(i-1)=b(i-1)/a(i-1);
    a(i)=a(i)-b(i-1)*c(i-1);
end

end

function [x] = LUsolve(a, b, c, d)
%     Solve a linear system LUx=b using backward substitution.
%     Input:
%     a, b, c: Output from LUdecomp
%     d: The RHS term of the linear system, N-element array
%     Output:
%     x: Solution, N-element array

x = zeros(size(d),'like',d);   % Initialize the solution array
N = length(d);
% ----------------------------------------
% TODO: Implement the backward substitution.
% ----------------------------------------
x(1)=d(1);
for j=2:N
    x(j)=d(j)-b(j-1)*x(j-1);
    x(N)=x(N)/a(N);
end
for k=N-1:-1:1
    x(k)=(x(k)-c(k)*x(k+1))/a(k);
end

end

function [y,t,U] = heatSolver(dy, dt, Q2)
%     Solves the unsteady heat transfer problem.
%     Input:
%     dy: Step size in y-direction, a float
%     dt: Time step size, a float
%     Q2: Heat flux of the cooling system, default Q2=5.0
%     Output:
%     y: Grid points in y-direction, (Ny+1)-element vector
%     t: All the time steps, (Nt+1)-element vector
%     U: An array containing all solutions, (_Nt+1)x(_Ny+1) matrix

% ----------------------------------------
% TODO: Comment on the lines below or develop your own implementation.
% ----------------------------------------

h  = 1.0;
Ti  = 10.0;

Ny = int16(ceil(h/dy));  % Determine the number of grid points
Nt = int16(ceil(Ti/dt));  % Determine the number of time steps
y  = linspace(0, h, Ny+1);  % Generate the grid point vector
t  = linspace(0, Ti, Nt+1);  % Generate the time step vector
U  = zeros(Nt+1, Ny+1);      % Allocate the array for numerical solutions

% Initialize the numerical discretization
[F, T, S, B1, B2] = initProblem(Ny, dy, dt, Q2);
% LU decomposition of the T matrix
[a, b, c] = LUdecomp(T);

for i = 1:Nt
    u = U(i,2:Ny)'; %Filling intermidiate matrix
    U(i+1,2:Ny) = LUsolve(a, b, c, S*u+F);
    U(i+1,1) = B1(1) + B1(2)*U(i+1,2) + B1(3)*U(i+1,3);
    U(i+1,Ny+1) = B2(1) + B2(2)*U(i+1, length(U(i+1,:))-1) + B2(3)*U(i+1, length(U(i+1,:))-2);
end

end

function [u] = anaSol(y, t, n, Q2)
%     Generates analytical solution
%     Input:
%     y: The grid points at which the analytical solutions are evaluated, N-element vector
%     t: The time at which the analytical solutions are evaluated, a float
%     n: Number of eigenfunctions to use, an integer
%     Q2: Heat flux of the cooling system
%     Output:
%     u: Analytical solutions evaluated at grid points *y* and time *t*, N-element vector
nu=1.0;
h=1;
Q1=10;
u =  zeros(size(y),'like',y);
% ----------------------------------------
% TODO: Implement the analytical solution
% ----------------------------------------
w= (Q1-Q2)/(2*h)*y.^2-Q1*y;
u=u+w;
u=u+(2*Q1+Q2)*h/6;
for i=1 :n+1
    p=i*pi/h;
    if mod(i,2)==0
        An=(2*h*(Q2-Q1))/(i*pi)^2;
    else
        An=(2*h*(-Q1-Q2))/(i*pi)^2;
    end
    u=u+An*cos(p*y)*exp(-nu*p^2*t);
end
u=u+nu*(Q1-Q2)/h*t;
end


