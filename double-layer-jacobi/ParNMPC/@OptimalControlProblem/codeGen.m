function codeGen(OCP)
    % init C if it is empty when the names of u x p have been changed
    if OCP.dim.mu == 0
%        OCP.C = symfun(zeros(0,1),[OCP.u;OCP.x;OCP.p]);
       OCP.C = sym(zeros(0,1));
    end
    % init G if it is empty
    if OCP.dim.z == 0
       OCP.G = symfun(zeros(0,1),[OCP.u;OCP.x;OCP.p]);
    end
    % 
    if OCP.isMEnabled
        % forced to 'Euler'
        OCP.discretizationMethod = 'Euler';
    else
        % init M if it is disabled
        OCP.M = sym(eye(OCP.dim.x));
    end
    showInfo(OCP);
    createNonemptySolution_FuncGen(OCP);
    %%
    disp('Generating OCP...');
    UXP = {OCP.u;OCP.x;OCP.p};
    %% LBarrier, LBarrieru, LBarrierx
%     LBarrier = symfun(OCP.LBarrier,[OCP.u;OCP.x;OCP.p]);
    LBarrier  = OCP.LBarrier;
    LBarrieru = jacobian(LBarrier,OCP.u);
    LBarrierx = jacobian(LBarrier,OCP.x);
    matlabFunction(LBarrier,LBarrieru,LBarrierx,...
        'File','./funcgen/OCP_GEN_LB_LBu_LBx',...
        'Vars',UXP,...
        'Outputs',{'LB','LBu','LBx'},'Optimize',false);
    matlabFunction(LBarrier,...
        'File','./funcgen/OCP_GEN_LB',...
        'Vars',UXP,...
        'Outputs',{'LB'},'Optimize',false);
    %% G
    matlabFunction(OCP.G,...
    'File','./funcgen/OCP_GEN_G',...
    'Vars',UXP,...
    'Outputs',{'G'},'Optimize',false);
    %% L, Lu, Lx
    Lu = jacobian(OCP.L,OCP.u);
    Lx = jacobian(OCP.L,OCP.x);
    matlabFunction(OCP.L,Lu,Lx,...
        'File','./funcgen/OCP_GEN_L_Lu_Lx',...
        'Vars',UXP,...
        'Outputs',{'L','Lu','Lx'},'Optimize',false);
    matlabFunction(OCP.L,...
    'File','./funcgen/OCP_GEN_L',...
    'Vars',UXP,...
    'Outputs',{'L'},'Optimize',false);
    %% C, Cu, Cx
    parIdx = sym('parIdx');
    UXPParIdx = {OCP.u;OCP.x;OCP.p;parIdx};

    if isa(OCP.C,'char')
        % external
        isExistCWrapper = exist('./C_Wrapper.m','file');
        if isExistCWrapper ~= 2
            copyfile('../ParNMPC/Wrapper/C_Wrapper.m');
            disp('Please specify your own C(u,x,p) function in C_Wrapper.m');
        else
            disp('C_Wrapper.m already exists and will be kept');
        end
        isExistC_Cu_Cx_Wrapper = exist('./C_Cu_Cx_Wrapper.m','file');
        if isExistC_Cu_Cx_Wrapper ~= 2
            copyfile('../ParNMPC/Wrapper/C_Cu_Cx_Wrapper.m');
            disp(['Please specify your own C_Cu_Cx(u,x,p) function in C_Cu_Cx_Wrapper.m', ...
                    ' (finite difference is used by default)']);   
        else
            disp('C_Cu_Cx_Wrapper.m already exists and will be kept');   
        end
        OCP.OCP_GEN_C_FuncGen();
        OCP.OCP_GEN_C_Cu_Cx_FuncGen();
    else
        if OCP.dim.mu == 0
            Cu = sym(zeros(0,OCP.dim.u));
            Cx = sym(zeros(0,OCP.dim.x));
        else
            Cu = jacobian(OCP.C,OCP.u);
            Cx = jacobian(OCP.C,OCP.x);
        end
        matlabFunction(OCP.C,Cu,Cx,...
        'File','./funcgen/OCP_GEN_C_Cu_Cx',...
        'Vars',UXPParIdx,...
        'Outputs',{'C','Cu','Cx'},'Optimize',false);
        matlabFunction(OCP.C,...
        'File','./funcgen/OCP_GEN_C',...
        'Vars',UXPParIdx,...
        'Outputs',{'C'},'Optimize',false);
    end

    %% f, fu, fx, M(optional)
    parIdx = sym('parIdx');
    UXPParIdx = {OCP.u;OCP.x;OCP.p;parIdx};

    if isa(OCP.f,'char')
        % external
        isExistfWrapper = exist('./f_Wrapper.m','file');
        if isExistfWrapper ~= 2
            copyfile('../ParNMPC/Wrapper/f_Wrapper.m');
            disp('Please specify your own f(u,x,p) function in f_Wrapper.m');
        else
            disp('f_Wrapper.m already exists and will be kept');
        end
        isExistf_fu_fx_Wrapper = exist('./f_fu_fx_Wrapper.m','file');
        if isExistf_fu_fx_Wrapper ~= 2
            copyfile('../ParNMPC/Wrapper/f_fu_fx_Wrapper.m');
            disp(['Please specify your own f_fu_fx(u,x,p) function in f_fu_fx_Wrapper.m', ...
                    ' (finite difference is used by default)']);   
        else
            disp('f_fu_fx_Wrapper.m already exists and will be kept');   
        end
        OCP.OCP_GEN_fdt_FuncGen();
        OCP.OCP_GEN_fdt_fudt_fxdt_FuncGen();
    else
        fdt  = OCP.f*OCP.deltaTau;
        fudt = jacobian(fdt,OCP.u);
        fxdt = jacobian(fdt,OCP.x);

        matlabFunction(fdt,fudt,fxdt,...
            'File','./funcgen/OCP_GEN_fdt_fudt_fxdt',...
            'Vars',UXPParIdx,...
            'Outputs',{'fdt','fudt','fxdt'},'Optimize',false);
        matlabFunction(fdt,...
            'File','./funcgen/OCP_GEN_fdt',...
            'Vars',UXPParIdx,...
            'Outputs',{'fdt'},'Optimize',false);
        
        % Generate Fu*v 
        vu = sym('vu',[OCP.dim.u,1]);
        Fu_v = fudt*vu;
        UXPvu = {OCP.u;OCP.x;OCP.p;vu};
        matlabFunction(Fu_v,...
            'File','./funcgen/OCP_GEN_Fuv',...
            'Vars',UXPvu,...
            'Outputs',{'Fu_v'},'Optimize',false);
        
        % Generate FuT*v
        vx = sym('vx',[OCP.dim.x,1]);
        FuT_v = fudt.'*vx;
        UXPvx = {OCP.u;OCP.x;OCP.p;vx};
        matlabFunction(FuT_v,...
            'File','./funcgen/OCP_GEN_FuTv',...
            'Vars',UXPvx,...
            'Outputs',{'FuT_v'},'Optimize',false);
        
        
        % Generate Fx_diag
        Fx = fxdt - eye(OCP.dim.x);
        Fx_diag = diag(Fx);
        matlabFunction(Fx_diag,...
            'File','./funcgen/OCP_GEN_Fx_diag',...
            'Vars',UXP,...
            'Outputs',{'Fx_diag'},'Optimize',false);
        
        % Generate (Fx - Fx_diag)*v
        FxmFx_diagvx = (Fx - diag(Fx_diag))*vx;
        matlabFunction(FxmFx_diagvx,...
            'File','./funcgen/OCP_GEN_FxmFx_diagvx',...
            'Vars',UXPvx,...
            'Outputs',{'FxmFx_diagvx'},'Optimize',false);
        
        % Generate (FxT - Fx_diag)*v
        FxTmFx_diagvx = (Fx.' - diag(Fx_diag))*vx;
        matlabFunction(FxTmFx_diagvx,...
            'File','./funcgen/OCP_GEN_FxTmFx_diagvx',...
            'Vars',UXPvx,...
            'Outputs',{'FxTmFx_diagvx'},'Optimize',false);
        
        
        % Generate FxT*v
        FxT_v = Fx.'*vx;
        matlabFunction(FxT_v,...
            'File','./funcgen/OCP_GEN_FxTv',...
            'Vars',UXPvx,...
            'Outputs',{'FxT_v'},'Optimize',false);

        
    end
    matlabFunction(OCP.M,...
        'File','./funcgen/OCP_GEN_M',...
        'Vars',UXP,...
    'Outputs',{'M'});
    if OCP.isMEnabled
        % Mx
        Mx = sym('Mx',[OCP.dim.x (OCP.dim.x*OCP.dim.x)]);
        for i = 1:OCP.dim.x
            Mx(:,(i-1)*OCP.dim.x+1:i*OCP.dim.x) = formula(diff(OCP.M,OCP.x(i)));
        end
        MxSymfun = symfun(Mx,[OCP.u;OCP.x;OCP.p]);
        matlabFunction(MxSymfun,...
            'File','./funcgen/OCP_GEN_Mx.m',...
            'Vars',UXP,...
            'Outputs',{'Mx'});
        % Mu
        Mu = sym('Mu',[OCP.dim.x (OCP.dim.x*OCP.dim.u)]);
        for i = 1:OCP.dim.u
            Mu(:,(i-1)*OCP.dim.x+1:i*OCP.dim.x) = formula(diff(OCP.M,OCP.u(i)));
        end
        MuSymfun = symfun(Mu,[OCP.u;OCP.x;OCP.p]);
        matlabFunction(MuSymfun,...
            'File','./funcgen/OCP_GEN_Mu.m',...
            'Vars',UXP,...
            'Outputs',{'Mu'});
    end
    disp('Done!');
end