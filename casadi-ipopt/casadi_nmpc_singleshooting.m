% addpath('C:\Users\dengh\Documents\casadi\casadi-windows-matlabR2016a-v3.5.5')
addpath('C:\Users\DengHaoyang\Documents\casadi-windows-matlabR2016a-v3.5.5')
clc
clear all
import casadi.*

M = 13;
numPeltier = 16;
nu = numPeltier;
nx = M*M - numPeltier;

Ts = 5; %[s]
N = 20; % prediction horizon

%% dynamics 
x = SX.sym('x', nx);
u = SX.sym('u', nu);

% Give names to x, u, p
T = SX.sym('T', M, M);
xCnt = 1;
uCnt = 1;
for yi=1:M
    for xi=1:M
        if mod(xi-1,4)==0 && mod(yi-1,4)==0
            T(xi,yi) = u(uCnt);
            uCnt = uCnt + 1;
        else
            T(xi,yi) = x(xCnt);
            xCnt = xCnt + 1;
        end
    end
end

% constants
k = 400; % thermal conductivity of copper, W/(m-K)
rho = 8960; % density of copper, kg/m^3
specificHeat = 386; % specific heat of copper, J/(kg-K)
thick = .01; % plate thickness in meters
emiss = .5; % emissivity of the plate surface
stefanBoltz = 5.670373e-8; % Stefan-Boltzmann constant, W/(m^2-K^4)
hCoeff = 1; % Convection coefficient, W/(m^2-K)
ta = 300; % The ambient temperature is assumed to be 300 degrees-Kelvin.
Cp = specificHeat;
tz = thick;
epsilon = emiss;
sigma = stefanBoltz;
hc = hCoeff;
Ta = ta;

width  = 1; 
height = 1;
h = width/M;

dTdt = SX.sym('dTdt', M,M);
for yi=1:M
    for xi=1:M
        if xi == 1
            TLeft   = T(xi,yi);
            TRright = T(xi+1,yi);
        elseif xi== M
            TLeft   = T(xi-1,yi);
            TRright = T(xi,yi);
        else
            TLeft   = T(xi-1,yi);
            TRright = T(xi+1,yi);
        end
        if yi == 1
            TUp   = T(xi,yi+1);
            TDown = T(xi,yi);
        elseif yi == M
            TUp   = T(xi,yi);
            TDown = T(xi,yi-1);
        else
            TUp   = T(xi,yi+1);
            TDown = T(xi,yi-1);
        end
        
        nabla2T = (TLeft + TRright + TUp + TDown - 4*T(xi,yi))/(h^2);
        
        Qc      = hc*(T(xi,yi) - Ta);
        Qr      = epsilon*sigma*(T(xi,yi)^4 - Ta^4);
        dTdt(xi,yi)    = (k*tz*nabla2T - 2*Qc - 2*Qr)/(rho*Cp*tz);
    end
end

fCnt = 1;
f = SX.sym('f', nx, 1);
for yi=1:M
    for xi=1:M
        if ~(mod(xi-1,4)==0 && mod(yi-1,4)==0)
          f(fCnt) =   dTdt(xi,yi);
          fCnt = fCnt + 1;
        end
    end
end
ffunc = Function('f',{x,u},{f}); % nonlinear mapping function f(x,u)
%% OCP
U = SX.sym('U', nu, N); 
X = SX.sym('X', nx, N+1);
P = SX.sym('P', M*M + nx); % parameters (which include the initial state and the reference)

obj = 0; % Objective function
Q = 1*eye(nx);
R = 0.1*eye(nu);
xRef = Ta + P(1:nx);
uRef = Ta + P(nx+1:nx+nu);

g = [];  % Eq constraints
% x0  = P(M*M + 1:end); % initial state
% g = [g; x0 - P(M*M + 1:end)]; 
X(:, 1) = P(M*M + 1:end);
for k = 1:N
    uk = U(:,k);
    % implicit euler based dynamics
%     rhs = X(:,k) + Ts * ffunc(xNext, uk);
    % explicit euler based dynamics
    X(:,k+1) = X(:,k) + Ts * ffunc(X(:,k), uk);  
%     g = [g; rhs - xNext];
    % objective
    L = 0.5*(X(:,k+1)-xRef).'*Q*(X(:,k+1)-xRef) + 0.5*(uk-uRef).'*R*(uk-uRef);
    obj = obj + L;
end
%%
TRefAllData   = load('TRefAll.mat');
TRefAll    = TRefAllData.TRefAll;
uRef = zeros(16,1);
TRef = zeros(153,1);
TRefCnt = 1;
uRefCnt = 1;
for yi=1:M
    for xi=1:M
        if ~(mod(xi-1,4)==0 && mod(yi-1,4)==0)
            TRef(TRefCnt) =   TRefAll(xi,yi);
            TRefCnt = TRefCnt + 1;
        else
            uRef(uRefCnt) =  TRefAll(xi,yi);
            uRefCnt = uRefCnt + 1;
        end
    end
end
%% solver
% X, U
OPT_variables = [reshape(U, nu*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', [], 'p', P);

opts = struct;
opts.ipopt.max_iter = 50;
opts.ipopt.print_level = 0; % 0, 3, 5
opts.ipopt.linear_solver = 'mumps';
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-3;
opts.ipopt.acceptable_obj_change_tol = 1e-6;
opts.ipopt.mu_target = 1e2;
opts.ipopt.mu_init = 1e2;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;
% Equality constraints
% args.lbg(1:nx*(N+1)) = 0;  
% args.ubg(1:nx*(N+1)) = 0;
% lower and upper
args.lbx(1:nu*N) = Ta;
args.ubx(1:nu*N)  = Ta + 400;
% parameters
x0 = (Ta)*ones(nx,1); % initial state
ref = [TRef; uRef]; % ref
args.p   = [ref;x0]; % set the values of the parameters vector
% initial guess
args.x0  = [repmat((Ta + 200) * ones(nu, 1), N,1)];
% solve
tic
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
             'lbg', [], 'ubg', [],'p', args.p);
toc
% get solution
uOptimal = reshape(full(sol.x(:)).', nu, N); 

plot(uOptimal.');
solver.stats.iter_count

%% closed-loop sim



