function xNext = SIM_Plant_RK4(u,x,p,Ts)

f = SIM_GEN_f(u,x,p);
k1 = Ts*(f);

f = SIM_GEN_f(u,x+k1/2,p);
k2 = Ts*(f);

f = SIM_GEN_f(u,x+k2/2,p);
k3 = Ts*(f);

f = SIM_GEN_f(u,x+k3,p);
k4 = Ts*(f);

xNext = x + (k1+2*k2+2*k3+k4)/6;