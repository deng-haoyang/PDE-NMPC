clear all
addpath('../ParNMPC/')
%% Formulate an OCP using Class OptimalControlProblem
M = 13;
numPeltier = 16;
% Create an OptimalControlProblem object
OCP = OptimalControlProblem(numPeltier,... % dim of inputs 
                            M*M - numPeltier,... % dim of states 
                            M*M + 1,... % dim of parameters
                            20);  % N: num of discritization grids

% Give names to x, u, p
T = sym(zeros(M,M));
xCnt = 1;
uCnt = 1;
for yi=1:M
    for xi=1:M
        if mod(xi-1,4)==0 && mod(yi-1,4)==0
            T(xi,yi) = OCP.u(uCnt);
            uCnt = uCnt + 1;
        else
            T(xi,yi) = OCP.x(xCnt);
            xCnt = xCnt + 1;
        end
    end
end

% Set the prediction horizon T
OCP.setT(OCP.p(end));

% Set the dynamic function f
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

dTdt = sym(zeros(M,M));
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
f = sym(zeros(OCP.dim.x,1));
for yi=1:M
    for xi=1:M
        if ~(mod(xi-1,4)==0 && mod(yi-1,4)==0)
          f(fCnt) =   dTdt(xi,yi);
          fCnt = fCnt + 1;
        end
    end
end

OCP.setf(f);
OCP.setDiscretizationMethod('Euler');

% Set the cost function L
Q = 1*eye(OCP.dim.x);%diag(OCP.p(1:OCP.dim.x));
R = 0.1*eye(OCP.dim.u);
xRef = Ta + OCP.p(1:OCP.dim.x);
uRef = Ta + OCP.p(OCP.dim.x+1:OCP.dim.x+OCP.dim.u);

L =    0.5*(OCP.x-xRef).'*Q*(OCP.x-xRef)...
     + 0.5*(OCP.u-uRef).'*R*(OCP.u-uRef);
OCP.setL(L);

% Set the linear constraints G(u,x,p)>=0
G =[OCP.u - Ta;...
   -OCP.u + Ta + 400];
OCP.setG(G);

% Generate necessary files
OCP.codeGen();
%% Configrate the solver using Class NMPCSolver

% Create a NMPCSolver object
nmpcSolver = NMPCSolver(OCP);

% Configurate the Hessian approximation method
nmpcSolver.setHessianApproximation('Newton');

% Generate necessary files
nmpcSolver.codeGen();
%% Solve the very first OCP for a given initial state and given parameters using Class OCPSolver

% Set the initial state
x0 =   (Ta)*ones(OCP.dim.x,1);

% Set the parameters
dim      = OCP.dim;
N        = OCP.N;
p      = zeros(dim.p,N);
p(1:end-1,:) = 0;   % setpoint
p(end,:)     = 100;    %T

% load TRef
TRefAllData   = load('TRefAll_gradient.mat');
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
for i=1:N
    p(1:153,i)   = TRef;
    p(154:169,i) = uRef;
end


% Solve the very first OCP 
solutionInitGuess.lambda = [randn(dim.lambda,1),zeros(dim.lambda,1)];
solutionInitGuess.mu     = randn(dim.mu,1);
solutionInitGuess.u      = (Ta+10)*ones(OCP.dim.u,1);
solutionInitGuess.x      = [x0,(Ta+10)*ones(OCP.dim.x,1)];
solutionInitGuess.z      = ones(dim.z,N);
solution = NMPC_SolveOffline(x0,p,solutionInitGuess,100,1000);

plot(solution.x([1 3 5],:).');

% Save to file
save GEN_initData.mat dim x0 p N

% Set initial guess
global ParNMPCGlobalVariable
ParNMPCGlobalVariable.solutionInitGuess = solution;
%% Define the controlled plant using Class DynamicSystem

% M(u,x,p) \dot(x) = f(u,x,p)
% Create a DynamicSystem object
plant = DynamicSystem(numPeltier,M*M - numPeltier,0);

% Set the dynamic function f
plant.setf(f); % same model 

% Generate necessary files
plant.codeGen();
