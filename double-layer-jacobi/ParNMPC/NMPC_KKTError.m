function [KKTError,costL,timeElapsed] = NMPC_KKTError(x0,p,rho,lambda,mu,u,x)

    timerStart = Timer();
    % Global variables
    global ParNMPCGlobalVariable
    discretizationMethod       = ParNMPCGlobalVariable.discretizationMethod;
    isMEnabled                 = ParNMPCGlobalVariable.isMEnabled;
    lambdaDim                  = ParNMPCGlobalVariable.dim.x;
    muDim                      = ParNMPCGlobalVariable.dim.mu;
    uDim                       = ParNMPCGlobalVariable.dim.u;
    xDim                       = ParNMPCGlobalVariable.dim.x;
    
    % init
    KKTError.stateEquation   = 0;
    KKTError.C               = 0;
    KKTError.Hu              = 0;
    KKTError.costateEquation = 0;
    
    [~,N]  = size(x);

    
    % Local
    KKTxEquation      = zeros(1, N);
    KKTC              = zeros(1, N);
    KKTHu             = zeros(1, N);
    KKTlambdaEquation = zeros(1, N);
    L                 = zeros(1,N);
        
%     parfor (i=1:1:DoP,numThreads)
    for i=1:1:N
        lambda_i = lambda(:,i);
        mu_i = mu(:,i);
        u_i = u(:,i);
        x_i = x(:,i);
        p_i = p(:,i);
                
        % Function and Jacobian
        [L_i,Lu_i, Lx_i]  = OCP_L_Lu_Lx(u_i,x_i,p_i);
        [~,LBu_i,LBx_i]   = OCP_LB_LBu_LBx(u_i,x_i,p_i);

%         [F_i,~,~] = OCP_F_Fu_Fx(u_i,x_i,p_i,discretizationMethod,isMEnabled,i);
        F_i = OCP_GEN_fdt(u_i,x_i,p_i) - x_i;
        % KKT
        if i > 1
            xPrev = x(:,i-1);
        else
            xPrev = x0;
        end
        if i < N
            lambdaNext = lambda(:,i+1);
        else
            lambdaNext = zeros(xDim,1);
        end
        xEq_i        = xPrev + F_i;
        % HuT_i      = Lu_i.' + Fu_i.'*lambda_i + rho*LBu_i.';
        HuT_i        = Lu_i.' + OCP_GEN_FuTv(u_i,x_i,p_i,lambda_i) + rho*LBu_i.';
        % HxT_i      = Lx_i.' + Fx_i.'*lambda_i + rho*LBx_i.';
        HxT_i        = Lx_i.' + OCP_GEN_FxTv(u_i,x_i,p_i,lambda_i) + rho*LBx_i.';
        lambdaEq_i  = lambdaNext + HxT_i;
        
        KKTxEquation(1,i) = norm(xEq_i,Inf);
        KKTC(1,i) = 0;
        KKTHu(1,i) = norm(HuT_i,Inf);
        KKTlambdaEquation(1,i) = norm(lambdaEq_i,Inf);

        L(:,i)               = L_i;
    end
    
    KKTError.stateEquation   = max(KKTxEquation(:));
    KKTError.C               = max(KKTC(:));
    KKTError.Hu              = max(KKTHu(:));
    KKTError.costateEquation = max(KKTlambdaEquation(:));
    costL = sum(L(:));
    timerEnd = Timer();
    timeElapsed = timerEnd  - timerStart;
end