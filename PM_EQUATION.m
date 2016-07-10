function [sol] = PM_EQUATION(v,M,gam)

% SOLVE PRANDTL-MEYER EQUATION
% Written by: JoshTheEngineer
% Started: 12/07/15
% Updated: 12/07/15 - Started code
%                   - Works as intended
% 
% PURPOSE
%   Solves the Prandtl-Meyer equation for either the Mach number or
%   the P-M angle, depending on the inputs
% 
% INPUTS
% - v   : Prandtl-Meyer angle [deg]
% - M   : Mach number []
% - gam : Ratio of specific heats []
% 
% OUTPUTS
% - sol : If M is input, then output is v [deg]
%         If v is input, then output is M []

% ---------------------- Solve for the P-M angle --------------------------
if (v == 0)
    term1 = sqrt((gam+1)/(gam-1));
    term2 = atand(sqrt((gam-1)*(M^2-1)/(gam+1)));
    term3 = atand(sqrt(M^2-1));
    
    sol = term1*term2 - term3;
end

% --------------------- Solve for the Mach number -------------------------
if (M == 0)
    problem.objective = @(M) (((sqrt(2.4/0.4))*(atand(sqrt((0.4)*(M^2-1)/(2.4)))))-...
                                    atand(sqrt(M^2-1)))-v;
    problem.x0        = [1 80];
    problem.solver    = 'fzero';
    problem.options   = optimset(@fzero);
    
    sol = fzero(problem);
end
