function [x,error] = jacobiIter(A,b,A_diag,numIter)
xTrue = A\b;
x      = b;
for i=1:numIter
    x = (b-(A-diag(A_diag))*x)./A_diag;
end

error = norm(x-xTrue);
end