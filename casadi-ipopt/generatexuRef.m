function xuRef =  generatexuRef(matfile)

uRef = zeros(16,1);
TRef = zeros(153,1);
TRefCnt = 1;
uRefCnt = 1;
% load TRef
TRefAllData   = coder.load(matfile);
TRefAll    = TRefAllData.TRefAll;
M = 13;

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

xuRef = [TRef;uRef];