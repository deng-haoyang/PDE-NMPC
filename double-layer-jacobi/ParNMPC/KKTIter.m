function [dlambda,du,dx] = KKTIter(Fx,Fu,Axx,Aux,Auu,KKT)

[xDim,uDim] = size(Fu);
KKTMatrix = [zeros(xDim,xDim),Fx,Fu;...
             Fx.',Axx,Aux.';...
             Fu.',Aux,Auu];
dAll = KKTMatrix\KKT;
dlambda =  dAll(1:xDim,1);
dx = dAll(xDim+1:2*xDim,1);
du = dAll(2*xDim+1:end,1);

end

