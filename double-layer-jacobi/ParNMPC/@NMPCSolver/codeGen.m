function codeGen(solver)
unknowns = [solver.OCP.lambda;...
            solver.OCP.mu;...
            solver.OCP.u;...
            solver.OCP.x;...
            solver.OCP.p];
        
rho = sym('rho',[1,1]);
LAll = solver.OCP.L + rho*solver.OCP.LBarrier;
% LAll = solver.OCP.L;

if  solver.OCP.isMEnabled == false
    switch solver.HessianApproximation
        case 'GaussNewton'
            A_ = LAll;
        case 'GaussNewtonLC'
            A_ = LAll + solver.OCP.mu.'*solver.OCP.C;
        case 'GaussNewtonLF'
            F = solver.OCP.f * solver.OCP.deltaTau - solver.OCP.x;
            A_ = LAll + solver.OCP.lambda.'*F;
        case 'Newton'
            F = solver.OCP.f * solver.OCP.deltaTau - solver.OCP.x;
            A_ = LAll + solver.OCP.mu.'*solver.OCP.C + solver.OCP.lambda.'*F;
        otherwise
            A_ = LAll;
    end
else
    switch solver.HessianApproximation
        case 'GaussNewton'
            A_ = LAll;
        case 'GaussNewtonLC'
            A_ = LAll + solver.OCP.mu.'*solver.OCP.C;
        otherwise
            A_ = LAll;
    end
end
showInfo(solver);
%% Generate Hessian for NMPC
disp('Generating Hessian...')
% A = symfun(A_,unknowns);
A = A_;
Au  = jacobian(A,solver.OCP.u);
Ax  = jacobian(A,solver.OCP.x);

Auu = jacobian(Au,solver.OCP.u);
Aux = jacobian(Au,solver.OCP.x);
Axx = jacobian(Ax,solver.OCP.x);
% % Condensing of z
% if solver.OCP.dim.z ~= 0
%     Gu = jacobian(solver.OCP.G,solver.OCP.u);
%     Gx = jacobian(solver.OCP.G,solver.OCP.x);
%     G  = solver.OCP.G;
%     AuuCondensed = formula(Auu) + formula(Gu.'*diag(solver.OCP.z./G)*Gu);%
%     AuxCondensed = formula(Aux) + formula(Gu.'*diag(solver.OCP.z./G)*Gx);%
%     AxxCondensed = formula(Axx) + formula(Gx.'*diag(solver.OCP.z./G)*Gx);%
% else
%     AuuCondensed = formula(Auu);%
%     AuxCondensed = formula(Aux);%
%     AxxCondensed = formula(Axx);%
% end
LambdaMuUXZP = {solver.OCP.lambda;...
               solver.OCP.mu;...
               solver.OCP.u;...
               solver.OCP.x;...
               solver.OCP.z;...
               solver.OCP.p;...
               rho};

% matlabFunction(AuuCondensed,AuxCondensed,AxxCondensed,...
%     'File','./funcgen/OCP_GEN_Auu_Aux_Axx_Condensed',...
%     'Vars',LambdaMuUXZP,...
%     'Outputs',{'AuuCondensed','AuxCondensed','AxxCondensed'});
matlabFunction(diag(Auu),...
    'File','./funcgen/OCP_GEN_Auu_diag',...
    'Vars',LambdaMuUXZP,...
    'Outputs',{'Auu_diag'},'Optimize',false);
matlabFunction(diag(Axx),...
    'File','./funcgen/OCP_GEN_Axx_diag',...
    'Vars',LambdaMuUXZP,...
    'Outputs',{'Axx_diag'},'Optimize',false);
disp('Done!');
end