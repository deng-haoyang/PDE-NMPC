function setL(OCP,L)
%     OCP.L = symfun(L,[OCP.u;OCP.x;OCP.p]);
    OCP.L = L;
end