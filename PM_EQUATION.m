% SOLVE PRANDTL-MEYER EQUATION
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/JoshTheEngineer
% Website   : www.joshtheengineer.com
% Started: 12/07/15
% Updated: 12/07/15 - Started code
%                   - Works as intended
%          11/??/17 - Updated to faster Mach number solver
% 
% PURPOSE
%   Solves the Prandtl-Meyer equation for either the Mach number or
%   the P-M angle, depending on the inputs
% 
% INPUTS
% - v : Prandtl-Meyer angle [deg]
% - M : Mach number []
% - g : Ratio of specific heats []
% 
% OUTPUTS
% - sol : If M is input, then output is v [deg]
%         If v is input, then output is M []

function [sol] = PM_EQUATION(v,M,g)

% For convenience
gm1 = g-1;
gp1 = g+1;

% ---------------------- Solve for the P-M angle --------------------------
if (v == 0)
    term1 = sqrt(gp1/gm1);
    term2 = atand(sqrt(gm1*(M^2-1)/gp1));
    term3 = atand(sqrt(M^2-1));
    
    sol = term1*term2 - term3;
end

% --------------------- Solve for the Mach number -------------------------
if (M == 0)
    dM  = 0.1;
    M   = 1;
    res = 1;
    while (res > 0.01)
        M2    = M + dM;
        funv1 = (-v*(pi/180)+(sqrt(gp1/gm1)*...
                    atan((sqrt(gm1*(M^2-1)/gp1)))-atan(sqrt(M^2-1))));
        funv2 = (-v*(pi/180)+(sqrt(gp1/gm1)*...
                    atan((sqrt(gm1*(M2^2-1)/gp1)))-atan(sqrt(M2^2-1))));
        dv_dm = (funv2-funv1)/dM;

        M   = M - funv1/dv_dm;
        res = abs(funv1);
    end
    sol = M;
end
