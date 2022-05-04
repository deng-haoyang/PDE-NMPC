function [dlambda,du,dx] = KKTIter(Fx,Fu,Axx,Aux,Auu,KKT)

[xDim,uDim] = size(Fu);
% KKTMatrix = [zeros(xDim,xDim),Fx,Fu;...
%              Fx.',Axx,Aux.';...
%              Fu.',Aux,Auu];
% 
% dAll = KKTMatrix\KKT;
% dlambda =  dAll(1:xDim,1);
% dx = dAll(xDim+1:2*xDim,1);
% du = dAll(2*xDim+1:end,1);

%%
Kx      = KKT(1:xDim,1);
Klambda = KKT(xDim+1:2*xDim,1);
Ku      = KKT(2*xDim+1:end,1);

[X1,X2] = FxAxxSolve(Fx,Axx,Kx,Klambda);
sys_du_b = Ku - [Fu.',Aux]*[X1;X2];

[X1,X2] = FxAxxSolve(Fx,Axx,Fu,Aux.');
sys_du_A = Auu - [Fu.',Aux]*[X1;X2];

% du = sys_du_A\sys_du_b;
[du,error] = jacobiIter(sys_du_A,sys_du_b,diag(Auu),2);


sys_dlambdax_b = [Kx;Klambda] - [Fu;Aux.']*du;
[dlambda,dx] = FxAxxSolve(Fx,Axx,sys_dlambdax_b(1:xDim),sys_dlambdax_b(xDim+1:end));

end

