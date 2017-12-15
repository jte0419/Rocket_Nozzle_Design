function [sol] = A_M_RELATION(ARatio,M,g,SubSup)
% =========================================================================
% - Solve for either the area ratio or Mach number
%   - Depends on what the user inputs are
% =========================================================================

% Set initial value of solution so it doesn't error out
sol = inf;

% Get and set convenient variables
gp1   = g + 1;
gm1   = g - 1;
gm1o2 = gm1/2;

% Solve for area ratio or Mach number
if (ARatio == 0)                                                            % SOLVE: AREA RATIO
    sol = sqrt((1/(M^2))*(((2/gp1)*(1+gm1o2*M^2)).^(gp1/gm1)));             % Solve
elseif (M == 0)                                                             % SOLVE: MACH NUMBER
    if (strcmpi(SubSup,'Sub'))                                              % Subsonic value
        problem.objective = @(M) sqrt((1/(M^2))*(((2/gp1)*...
                                    (1+gm1o2*M^2))^(gp1/gm1))) - ARatio;    % Objective function
        problem.x0        = [1e-6 1];                                       % Solver bounds
        problem.solver    = 'fzero';                                        % Find the zero
        problem.options   = optimset(@fzero);                               % Default options
        sol               = fzero(problem);                                 % Solve
    elseif (strcmpi(SubSup,'Sup'))                                          % Supersonic value
        problem.objective = @(M) sqrt((1/(M^2))*(((2/gp1)*...
                                    (1+gm1o2*M^2))^(gp1/gm1))) - ARatio;    % Objective function
        problem.x0        = [1 50];                                         % Solver Bounds
        problem.solver    = 'fzero';                                        % Find the zero
        problem.options   = optimset(@fzero);                               % Default options
        sol               = fzero(problem);                                 % Solve
    end
end

if (sol == inf)
    fprintf('Error in A_M_RELATION function\n\n');
end