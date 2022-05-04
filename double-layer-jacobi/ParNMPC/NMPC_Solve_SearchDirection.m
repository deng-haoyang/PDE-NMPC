function [lambda,mu,u,x,z,KKTError,costL,timeElapsed] = ...
                NMPC_Solve_SearchDirection(x0,p,rho,lambda,mu,u,x,z) %#codegen
    tStart = Timer();
    
    % global variables
    global ParNMPCGlobalVariable
    discretizationMethod       = ParNMPCGlobalVariable.discretizationMethod;
    isMEnabled                 = ParNMPCGlobalVariable.isMEnabled;
    nonsingularRegularization  = ParNMPCGlobalVariable.nonsingularRegularization;
    descentRegularization      = ParNMPCGlobalVariable.descentRegularization;
    isApproximateInvFx         = ParNMPCGlobalVariable.isApproximateInvFx;
    lambdaDim                  = ParNMPCGlobalVariable.dim.x;
    muDim                      = ParNMPCGlobalVariable.dim.mu;
    uDim                       = ParNMPCGlobalVariable.dim.u;
    xDim                       = ParNMPCGlobalVariable.dim.x;
    zDim                       = ParNMPCGlobalVariable.dim.z;
    
    descentRegularization = 0.5;
    
    % parallel seg
    [~,sizeSeg,DoP]  = size(u);
    [~,N]  = size(u);

    if coder.target('MATLAB') % Normal excution
        % in serial
        numThreads = 0;
    else % Code generation
        numThreads = DoP;
    end
    
    % backup 
    lambda_k = lambda;
    mu_k     = mu;
    u_k      = u;
    x_k      = x;
    z_k      = z;
    
    % line search parameters
    stepSizeZ     = ones(DoP,1);
    stepSizeG     = ones(DoP,1);
    
    % regularization
%     NSRMatrix  = -nonsingularRegularization*eye(muDim);
%     uDRMatrix  =  descentRegularization*eye(uDim);
%     xDRMatrix  =  descentRegularization*eye(xDim);
    
    % local variables
    lambdaNext      = zeros(lambdaDim,sizeSeg,DoP);
    xPrev           = zeros(xDim,sizeSeg,DoP);

%     dlambda         = zeros(lambdaDim,sizeSeg,DoP);
%     dx              = zeros(xDim,sizeSeg,DoP);
    
%     LAMBDA          = zeros(xDim,xDim,sizeSeg,DoP);
    KKTxEquation      = zeros(1, DoP);
    KKTC              = zeros(1, DoP);
    KKTHu             = zeros(1, DoP);
    KKTlambdaEquation = zeros(1, DoP);
    KKTError.stateEquation   = 0;
    KKTError.C               = 0;
    KKTError.Hu              = 0;
    KKTError.costateEquation = 0;
    
    % coupling variable for each segment
    xPrev(:,1,1) = x0;
    for i=2:1:DoP
        xPrev(:,1,i) = x(:,sizeSeg,i-1);
        lambdaNext(:,sizeSeg,i-1) = lambda(:,1,i);
    end
    
    lambdaCorrected = zeros(xDim,N);

    F        = zeros(xDim,N);
    HxT      = zeros(xDim,N);
    xEq      = zeros(xDim,N);
    C        = zeros(muDim,N);
    HuT      = zeros(uDim,N);
    lambdaEq = zeros(lambdaDim,N);
    L        = zeros(1,N);
    
    Auu_diag = zeros(uDim,N);
    Axx_diag = zeros(xDim,N);
    Fx_diag  = zeros(xDim,N);
    %% Backward
    for i=N:-1:1
        lambda_i = lambda(:,i);
        mu_i     = mu(:,i);
        u_i      = u(:,i);
        x_i = x(:,i);
        z_i = z(:,i);
        p_i = p(:,i);
        
        % Function and Jacobian
        [L_i,Lu_i, Lx_i]  = OCP_L_Lu_Lx(u_i,x_i,p_i);
        [~,LBu_i,LBx_i]   = OCP_LB_LBu_LBx(u_i,x_i,p_i);

%         [F_i,~,~] = OCP_F_Fu_Fx(u_i,x_i,p_i,discretizationMethod,isMEnabled,i);
        F_i = OCP_GEN_fdt(u_i,x_i,p_i) - x_i;
%         % Hessian
%         % --- AuuCondensed_j_i = Auu_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gu_j_i);
%         % --- AuxCondensed_j_i = Aux_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
%         % --- AxxCondensed_j_i = Axx_j_i + Gx_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
%         [Auu_i,Aux_i,Axx_i] = OCP_Auu_Aux_Axx_Condensed(lambda_i,mu_i,u_i,x_i,z_i,p_i);
%         Axx_i   = Axx_i + xDRMatrix; % descent regularization
%         Auu_i   = Auu_i + uDRMatrix; % descent regularization
        
        Fx_diag_i  = OCP_GEN_Fx_diag(u_i,x_i,p_i);
        Axx_diag_i = OCP_GEN_Axx_diag(lambda_i,mu_i,u_i,x_i,z_i,p_i,rho) + 0;
        Auu_diag_i = OCP_GEN_Auu_diag(lambda_i,mu_i,u_i,x_i,z_i,p_i,rho) + descentRegularization;
        
        % KKT
        if i > 1
            xPrev = x(:,i-1);
        else
            xPrev = x0;
        end
        if i < N
            lambdaNext = lambdaCorrected(:,i+1);
        else
            lambdaNext = zeros(xDim,1);
        end
        xEq_i      = xPrev + F_i;
        % HuT_i      = Lu_i.' + Fu_i.'*lambda_i + rho*LBu_i.';
        HuT_i        = Lu_i.' + OCP_GEN_FuTv(u_i,x_i,p_i,lambda_i) + rho*LBu_i.';
        % HxT_i      = Lx_i.' + Fx_i.'*lambda_i + rho*LBx_i.';
        HxT_i        = Lx_i.' + OCP_GEN_FxTv(u_i,x_i,p_i,lambda_i) + rho*LBx_i.';
        lambdaEq_i  = lambdaNext + HxT_i;
        
        % iteration 
        KKT_i = [xEq_i;lambdaEq_i;HuT_i];
        
        [dlambda_i,~,~] = KKTIter_2dheat(Fx_diag_i,Axx_diag_i,Auu_diag_i,lambda_i,mu_i,u_i,x_i,z_i,p_i,KKT_i);
%         [dlambda_i,~,~] = KKTIter(Fx_i,Fu_i,Axx_i,Aux_i,Auu_i,KKT_i);
        
        lambdaCorrected(:,i) = lambda(:,i) -  dlambda_i;
        
        % save
        F(:,i)   = F_i;
        HxT(:,i) = HxT_i;
        HuT(:,i) = HuT_i;
        L(:,i)   = L_i;
        xEq(:,i)   = xEq_i;
        lambdaEq(:,i)   = lambdaEq_i;
        
        Auu_diag(:,i) = Auu_diag_i;
        Axx_diag(:,i) = Axx_diag_i;
        Fx_diag(:,i)  = Fx_diag_i;
    end
    %% Forward
    for i=1:1:N
        % iteration 
        F_i      =  F(:,i);
        HxT_i    =  HxT(:,i);
        HuT_i    =  HuT(:,i);
        
        Auu_diag_i =  Auu_diag(:,i);
        Axx_diag_i =  Axx_diag(:,i);
        Fx_diag_i  =  Fx_diag(:,i);
        
        if i > 1
            xPrev = x(:,i-1);
        else
            xPrev = x0;
        end
        if i < N
            lambdaNext = lambdaCorrected(:,i+1);
        else
            lambdaNext = zeros(xDim,1);
        end
        xEq_i      = xPrev + F_i;
        lambdaEq_i  = lambdaNext + HxT_i;
        
        % iteration
        KKT_i = [xEq_i;lambdaEq_i;HuT_i];
        
%         [dlambda_i,du_i,dx_i] = KKTIter(Fx_i,Fu_i,Axx_i,Aux_i,Auu_i,KKT_i);
        [dlambda_i,du_i,dx_i] = KKTIter_2dheat(Fx_diag_i,Axx_diag_i,Auu_diag_i,lambda_i,mu_i,u_i,x_i,z_i,p_i,KKT_i);

        lambda(:,i) = lambda(:,i) -  dlambda_i;
        x(:,i) = x(:,i) -  dx_i;
        u(:,i) = u(:,i) - du_i;
    end
    %% line search
    if zDim ~= 0
        % G feasibility
        G   = OCP_G(u,x,p);
        G_k = OCP_G(u_k,x_k,p);
        dG  = G - G_k;
        stepSizeG = (0.005 - 1).*(G_k./dG);
        stepSizeG(stepSizeG>1 | stepSizeG<0) = 1;
        stepSizeMaxG = min(stepSizeG(:));
        if isempty(stepSizeMaxG)
            stepSizeMaxG = 1;
        end
    end
    %
    stepSizePrimal = stepSizeMaxG;
    if stepSizePrimal ~= 1
        lambda = (1-stepSizePrimal)*lambda_k + stepSizePrimal* lambda;
        mu     = (1-stepSizePrimal)*mu_k     + stepSizePrimal* mu;
        u      = (1-stepSizePrimal)*u_k      + stepSizePrimal* u;
        x      = (1-stepSizePrimal)*x_k      + stepSizePrimal* x;
    end
%     stepSizeDual = stepSizeMaxZ;
%     if stepSizeDual ~= 1
%         z = (1-stepSizeDual)*z_k + stepSizeDual*z;
%     end

    %%
    KKTError.stateEquation   = norm(xEq(:),Inf);
    KKTError.C               = 0;
    KKTError.Hu              = norm(HuT(:),Inf);
    KKTError.costateEquation = norm(lambdaEq(:),Inf);
    costL = sum(L(:));
    %%
    tEnd = Timer();
    timeElapsed = tEnd  - tStart;
end
