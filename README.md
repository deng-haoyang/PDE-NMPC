# Requirements
- MATLAB Coder (with a supported compiler)
- Symbolic Math Toolbox

# Run double-lay-jacobi in MATLAB
1. Go to double-lay-jacobi/2dHeat/
2. Run NMPC_Problem_Formulation.m
3. Run Simu_Matlab_Codegen.m to generate C code
4. Run ./Simu_Matlab.exe in MATLAB Command Line

or 

3. Simu_Matlab.m
4. Run visualRec.m to visualize the simulation results


# Run ParNMPC in MATLAB
1. Go to ParNMPC/2dHeat/
2. Run NMPC_Problem_Formulation.m
3. Run Simu_Matlab_Codegen.m to generate C code
4. Run ./Simu_Matlab.exe in MATLAB Command Line

or 

3. Run Simu_Matlab.m
4. Run visualRec.m to visualize the simulation results


# Run casadi-ipopt in MATLAB
## Run multiple-shooting NMPC
1. Go to casadi-ipopt/2dHeat/
2. Run casadi_nmpc_multipleshooting.m to generate C code
3. Run Simu_Matlab_casadi_nmpc_multipleshooting.m

## Run single-shooting NMPC
1. Go to casadi-ipopt/2dHeat/
2. Run casadi_nmpc_singleshooting.m
3. Run Simu_Matlab_casadi_nmpc_singleshooting.m
