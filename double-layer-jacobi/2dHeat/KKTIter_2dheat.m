function [dlambda,du,dx] = KKTIter_2dheat(Fx_diag,Axx_diag,Auu_diag,lambda_i,mu_i,u_i,x_i,z_i,p_i,KKT)

xDim = length(x_i);

Kx      = KKT(1:xDim,1);
Klambda = KKT(xDim+1:2*xDim,1);
Ku      = KKT(2*xDim+1:end,1);
%
numIter_jacobi_Auu = 2;


du = Ku./Auu_diag;
for i=1:numIter_jacobi_Auu
    Fudu = OCP_GEN_Fuv(u_i,x_i,p_i,du);
    % [Fu;Axu]*du - [Kx;Klambda] (Axu = 0)
    sys_b = [Fudu;zeros(xDim,1)] - [Kx;Klambda];
    sys_x = FxAxxSolve_2dheat(Fx_diag,Axx_diag,lambda_i,mu_i,u_i,x_i,z_i,p_i,sys_b(1:xDim),sys_b(xDim+1:end));
    du = (Ku + OCP_GEN_FuTv(u_i,x_i,p_i,sys_x(1:xDim)))./Auu_diag;
end
% recover dlambda, dx
%Fu*du
Fudu = OCP_GEN_Fuv(u_i,x_i,p_i,du);
[dlambda,dx] = FxAxxSolve_2dheat(Fx_diag,Axx_diag,lambda_i,mu_i,u_i,x_i,z_i,p_i,Kx - Fudu,Klambda);
end
%%
function [x1,x2] = FxAxxSolve_2dheat(Fx_diag,Axx_diag,lambda_i,mu_i,u_i,x_i,z_i,p_i,b1,b2)
numIter_jacobi_Fx = 2;

% initial guess
x2 = b1./Fx_diag;
x1 = (b2 - Axx_diag.*x2)./Fx_diag;

for i=1:numIter_jacobi_Fx
    % (Fx-Fx_diag)*x2
    FxmFx_diagx2 = OCP_GEN_FxmFx_diagvx(u_i,x_i,p_i,x2);
    x2 = (b1 - FxmFx_diagx2)./Fx_diag;
end

for i=1:numIter_jacobi_Fx
    % (Fx.'-Fx_diag)*x1
    FxTmFx_diagx2 = OCP_GEN_FxTmFx_diagvx(u_i,x_i,p_i,x1);
    x1 = (b2 - FxTmFx_diagx2 - Axx_diag.*x2)./Fx_diag;
end

end

