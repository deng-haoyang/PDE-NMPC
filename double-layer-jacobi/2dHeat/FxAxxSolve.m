function [X1,X2] = FxAxxSolve(Fx,Axx,B1,B2)
    [xDim,~] = size(Fx);
    X2 = FxSolve(Fx,B1);
    X1 = FxTSolve(Fx.',B2 - Axx*X2);
end

function X = FxSolve(Fx,B)
%     X = Fx\B;
    [X,error] = jacobiIter(Fx,B,diag(Fx),10);
end


function X = FxTSolve(FxT,B)
%     X = FxT\B;
    [X,error] = jacobiIter(FxT,B,diag(FxT),10);
end