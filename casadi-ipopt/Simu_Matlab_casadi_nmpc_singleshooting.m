%% For closed-loop simulation or code generation
simuLength = 1000;   
Ts         = 5; % sampling interval
simuSteps  = floor(simuLength/Ts);

% Load data
x0 = (Ta)*ones(nx,1); % initial state

% define record variables
rec.x       = zeros(simuSteps+1,nx);
rec.x(1,:)  = x0.';
rec.u       = zeros(simuSteps,  nu);
rec.numIter = zeros(simuSteps,1);
rec.error   = zeros(simuSteps,1);
rec.cost    = zeros(simuSteps,1);
rec.t       = zeros(simuSteps,1);
rec.cpuTime                = zeros(simuSteps,1);
rec.cpuTimeSearchDirection = zeros(simuSteps,1);
rec.cpuTimeLineSearch      = zeros(simuSteps,1);
rec.cpuTimeKKTError        = zeros(simuSteps,1);
%% ref
xuRef_gradient = generatexuRef("TRefAll_gradient.mat");
xuRef_v2 = generatexuRef("TRefAll_v2.mat");
% initial guess
args.x0  = [repmat((Ta + 200) * ones(nu, 1), N,1)];
xuRef = xuRef_gradient;

%% initial guess
% args.p   = [xuRef;x0Measured]; % set the values of the parameters vector
% % solve
% sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
%              'lbg', args.lbg, 'ubg', args.ubg,'p', args.p);
% args.x0  = sol.x;
%% Simulation
for step = 1:simuSteps %simulation steps
    % Solve the optimal control problem
    x0Measured = x0;    
    % parameters
    args.p   = [xuRef;x0Measured]; % set the values of the parameters vector
    % solve
    tStart = tic;
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
                 'lbg', [], 'ubg', [], 'p', args.p);
    cpuTime = toc(tStart);
    uOptimalTraj = reshape(full(sol.x(:)).', nu, N); 
%     xOptimaltraj = reshape(full(sol.x(1:nx*(N+1))).', nx, N + 1); 
    
    % initial guess
    args.x0  = sol.x;

%     tTotal        = tTotal + output.timeElapsed.total;
    
    % Obtain the first optimal control input
    uOptimal = uOptimalTraj(:,1);
    % System simulation by the 4th-order Explicit Runge-Kutta Method
    pSimVal = zeros(0,1);
    x0 = SIM_Plant_RK4(uOptimal(:,1), x0, pSimVal, Ts);

    % Update parameters
    if step >= 101 && step<111
        ratio = (step-100)/10;
        xuRef = ratio*xuRef_v2 + (1-ratio)*xuRef;
    end
    
    % Record data
    rec.x(step+1,:)      = x0Measured.';
    rec.u(step,:)        = uOptimal.';
    rec.cpuTime(step,:)  = cpuTime;
    rec.t(step,:)        = step*Ts;
    rec.numIter(step,:)  = solver.stats.iter_count;
    if coder.target('MATLAB')
         disp(['Step: ',num2str(step),'/',num2str(simuSteps),...
               '   iterTotal: ', num2str(rec.numIter(step,:))])
    end
end
%% Log to file
if coder.target('MATLAB')% Normal excution
    save('GEN_log_rec.mat','rec');
end