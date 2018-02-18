% METHOD OF CHARACTERISTICS - NOZZLE DESIGN
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/JoshTheEngineer
% Website   : www.joshtheengineer.com
% Started: 06/24/16
% Updated: 06/24/16 - Started code
%                   - Works as expected, but no plotting yet
%          06/25/16 - Plotting works now
%          06/27/16 - Added comments
%          11/18/17 - Updated inputs
%                   - Fixed gamma issue
%                   - Outputs GMSH file (still need to mesh and export)

% PURPOSE
% - Solve for the isentropic nozzle geometry from the givens
% - Can either solve using defined exit Mach number or defined area ratio
% 
% GIVENS
% - g     : Specific heat ratio []
% - Ae_At : Area ratio of nozzle []
% - Me    : Exit Mach number []

clear;
clc;

%% INPUTS

% Geometry
Dstar      = 2;                                                             % Diameter [length]
geom       = 1;                                                             % [1] Point, [2] Circle
radius     = 1;                                                             % Radius of throat [length]
M1         = 1.0001;                                                        % Throat Mach number (> 1) []
numChar    = 50;                                                            % Number of characteristic lines to use [#]
g          = 1.4;                                                           % Specific heat ratio []
R          = 287;                                                           % Specific gas constant [J/kg*K]
outputMesh = 0;                                                             % Output mesh [1 = Yes, 0 = No]

% Set knowns
Me_Set    = 0;                                                              % Exit Mach number []
Ae_At_Set = 3;                                                              % Nozzle area ratio []

if (Me_Set == 0)
    Me_Set = A_M_RELATION(Ae_At_Set,0,g,'Sup');
elseif (Ae_At_Set == 0)
    Ae_At_Set = A_M_RELATION(0,Me_Set,g,'Sup');
end

% Get anle array based on maximum throat turn angle
thetaMax = PM_EQUATION(0,Me_Set,g)/2;                                       % Maximum throat expansion angle [deg]
theta    = linspace(0,thetaMax,numChar)';                                   % Throat turn angle array [deg]

% Nozzle specific properties
gm1o2 = (g-1)/2;
togp1 = 2/(g+1);
gogm1 = g/(g-1);

P0 = 7e6;
T0 = 3558;
Ps = P0*(togp1^gogm1);
Ts = T0*(togp1);
as = sqrt(g*R*Ts);
Pe = P0/((1+gm1o2*Me_Set^2)^gogm1);
Te = T0/(1+gm1o2*Me_Set^2);

fprintf('Ae/At = %2.3f \n',Ae_At_Set);
fprintf('Me    = %2.3f \n',Me_Set);
fprintf('P0 = %2.3f [Pa]\n',P0);
fprintf('P* = %2.3f [Pa]\n',Ps);
fprintf('Pe = %2.3f [Pa]\n',Pe);
fprintf('T0 = %2.3f [K]\n',T0);
fprintf('T* = %2.3f [K]\n',Ts);
fprintf('Te = %2.3f [K]\n',Te);
fprintf('a* = %2.3f [m/s]\n',as);

%% INITIAL SETUP OF KNOWNS

% Initialize solution variables
dataENS  = cell(numChar,numChar+1);                                         % Expansion and non-simple regions
dataS    = cell(numChar,1);                                                 % Straightening region

% ====================
% ==== EXPANSION =====
% ====================
for R = 1:1:numChar                                                         % Loop through all (-) characteristics
    dataENS{R,1}.M = 0;                                                     % Set all Mach numbers to zero
    if (R == 1)                                                             % For the first point
        dataENS{R,1}.M = M1;                                                % Set the Mach number from the throat M
    end
    dataENS{R,1}.theta = theta(R);                                          % Angle w.r.t horizontal [deg]
    dataENS{R,1}.nu    = 0;                                                 % Prandtl-Meyer angle [deg]
    dataENS{R,1}.mu    = 0;                                                 % Mach angle [deg]
    dataENS{R,1}.Kp    = 0;                                                 % Plus characteristic constant [deg]
    dataENS{R,1}.Km    = 0;                                                 % Minus characteristic constant [deg]
    dataENS{R,1}.X     = 0;                                                 % Point X location
    dataENS{R,1}.Y     = 0;                                                 % Point Y location
    dataENS{R,1}.dydx  = 0;                                                 % Point slope
    dataENS{R,1}.tsip  = 0;                                                 % Convenient avg parameter for (+) characteristic
    dataENS{R,1}.tsim  = 0;                                                 % Convenient avg parameter for (-) characteristic
    
    % Apply geometry
    if (geom == 1)                                                          % Zero-radius corner (just a point)
        dataENS{R,1}.X = 0;                                                 % All X-values are at the throat
        dataENS{R,1}.Y = Dstar/2;                                           % All Y-values are at half the throat diameter
    elseif (geom == 2)                                                      % Finite-radius corner
        dataENS{1,1}.X = 0;                                                 % First X-value is at the throat
        dataENS{1,1}.Y = Dstar/2;                                           % First Y-value is at half the throat diameter
        
        for i = 2:1:numChar                                                 % For the rest of the starting line characteristics
            dx = radius*sind(theta(i));                                     % Change in X depends on radius and angle
            dy = radius - radius*cosd(theta(i));                            % Change in Y depends on radius and angle
            
            dataENS{i,1}.X = dataENS{1,1}.X + dx;                           % Add the X change to the first X point
            dataENS{i,1}.Y = dataENS{1,1}.Y + dy;                           % Add the Y change to the first Y point
        end
    end
end

% =====================
% ==== NON-SIMPLE =====
% =====================
for R = 1:1:numChar                                                         % Loop over all (-) characteristics
    for L = 2:1:numChar+1                                                   % Loop over all (+) characteristics
        dataENS{R,L}.M     = 0;                                             % Mach number
        dataENS{R,L}.theta = 0;                                             % Flow angle w.r.t horizontal [deg]
        dataENS{R,L}.nu    = 0;                                             % Prandtl-Meyer angle [deg]
        dataENS{R,L}.mu    = 0;                                             % Mach angle [deg]
        dataENS{R,L}.Kp    = 0;                                             % Plus characteristic constant [deg]
        dataENS{R,L}.Km    = 0;                                             % Minus characteristic constant [deg]
        dataENS{R,L}.X     = 0;                                             % Point X location
        dataENS{R,L}.Y     = 0;                                             % Point Y location
        dataENS{R,L}.dydx  = 0;                                             % Point slope
        dataENS{R,L}.tsip  = 0;                                             % Convenient avg parameter for (+) characteristic
        dataENS{R,L}.tsim  = 0;                                             % Convenient avg parameter for (-) characteristic
    end
end

% ========================
% ==== STRAIGHTENING =====
% ========================
for L = 1:1:numChar                                                         % Loop through all (+) characteristics
    dataS{L}.M     = 0;                                                     % Mach number []
    dataS{L}.theta = 0;                                                     % Flow angle w.r.t. horizontal[deg]
    dataS{L}.nu    = 0;                                                     % Prandtl-Meyer angle [deg]
    dataS{L}.mu    = 0;                                                     % Mach angle [deg]
    dataS{L}.Kp    = 0;                                                     % Plus characteristic constant [deg]
    dataS{L}.Km    = 0;                                                     % Minus characteristic constant [deg]
    dataS{L}.X     = 0;                                                     % Point X location
    dataS{L}.Y     = 0;                                                     % Point Y location
    dataS{L}.dydx  = 0;                                                     % Point slope
    dataS{L}.tsi   = 0;                                                     % Convenient nozzle contour parameter
    dataS{L}.tsip  = 0;                                                     % Convenient avg parameter for (+) characteristic
    dataS{L}.tsim  = 0;                                                     % Convenient avg parameter for (-) characteristic
end

%% EXPANSION REGION

% Indicate that we are solving expansion region, start section timer
fprintf('Solving expansion region: ');                                      % Print to the command window
tic;                                                                        % Start the timer

% Starting point Prandtl-Meyer angle [deg]
dataENS{1,1}.nu = PM_EQUATION(0,dataENS{1,1}.M,g);                          % PM angle from the throat Mach number [deg]

% Loop through all characteristics (Rows)
for R = 1:1:numChar                                                         % Loop over all the (-) characteristics
    % All throat K+ characteristic constants
    dataENS{R,1}.Kp = dataENS{1,1}.theta - dataENS{1,1}.nu;                 % All the same - coming from throat [deg]
    
    % PM angles, Mach numbers, Mach angles, and K- characterstic constants
    if (R ~= 1)                                                             % If we are not on the first characteristic (values already defined)
        dataENS{R,1}.nu = dataENS{R,1}.theta - dataENS{R,1}.Kp;             % Prandtl-Meyer angle (don't need first characteristic)
        dataENS{R,1}.M  = PM_EQUATION(dataENS{R,1}.nu,0,g);                 % Mach number (don't need first characteristic)
    end
    dataENS{R,1}.mu = asind(1/dataENS{R,1}.M);                              % Mach angle [deg]
    dataENS{R,1}.Km = dataENS{R,1}.theta + dataENS{R,1}.nu;                 % Minus characteristic constant [deg]
end

% Print the time taken to solve the expansion region
fprintf('%2.2f [sec]\n',toc);                                               % Print the time

%% NON-SIMPLE REGION

% Indicate that we are solving non-simple region, start section timer
fprintf('Solving non-simple region: ');                                     % Print to the command window
tic;                                                                        % Start the timer

startR = 1;                                                                 % Index for starting value of (-) characteristic
for L = 2:1:numChar+1                                                       % Loop through all (+) characteristics
    for R = startR:1:numChar                                                % Loop through appropriate (-) characteristics
        
        if (L-1 == R)                                                       % If we are on the first (+) characteristic
            dataENS{R,L}.nu = dataENS{R,1}.Km;                              % Prandtl-Meyer angle from (-) constant [deg]
            dataENS{R,L}.M  = PM_EQUATION(dataENS{R,L}.nu,0,g);             % Mach number from PM equation using nu as input []
            dataENS{R,L}.mu = asind(1/dataENS{R,L}.M);                      % Mach angle [deg]
            dataENS{R,L}.Kp = dataENS{R,L}.theta - dataENS{R,L}.nu;         % Plus characteristic constant [deg]
            dataENS{R,L}.Km = dataENS{R,L}.theta + dataENS{R,L}.nu;         % Minus characteristic constant [deg]
        else                                                                % For all other (+) characteristics
            dataENS{R,L}.Kp    = dataENS{R-1,L}.Kp;                         % Plus characteristic constant [deg]
            dataENS{R,L}.Km    = dataENS{R,L-1}.Km;                         % Minus characteristic constant [deg]
            dataENS{R,L}.theta = 0.5*(dataENS{R,L}.Km + dataENS{R,L}.Kp);   % Angle w.r.t. horizontal [deg]
            dataENS{R,L}.nu    = 0.5*(dataENS{R,L}.Km - dataENS{R,L}.Kp);   % Prandtl-Meyer angle [deg]
            dataENS{R,L}.M     = PM_EQUATION(dataENS{R,L}.nu,0,g);          % Mach number []
            dataENS{R,L}.mu    = asind(1/dataENS{R,L}.M);                   % Mach angle [deg]
        end        
    end
    
    startR = startR + 1;                                                    % Increment the starting (-) characteristic counter
end

% Print the time taken to solve the non-simple region
fprintf('%2.2f [sec]\n',toc);                                               % Print the time

%% STRAIGHTENING REGION

% Indicate that we are solving straightening region, start section timer
fprintf('Solving straightening region: ');                                  % Print to the command window
tic;                                                                        % Start the timer

% Solve for the fully straightened variables
dataS{end}.theta = 0;                                                       % Flow is back to horizontal [deg]
dataS{end}.M     = dataENS{end,end}.M;                                      % Exit Mach number []
dataS{end}.nu    = PM_EQUATION(0,dataENS{end}.M,g);                         % Exit Prandtl-Meyer angle [deg]
dataS{end}.mu    = asind(1/dataS{end}.M);                                   % Exit Mach angle [deg]
dataS{end}.Kp    = dataS{end}.theta - dataS{end}.nu;                        % Exit plus characteristic constant [deg]
dataS{end}.Km    = dataS{end}.theta + dataS{end}.nu;                        % Exit minus characteristic constant [deg]

% Using known exit values, solve for the rest of the straightening region
for L = 1:1:numChar-1
    dataS{L}.Km    = dataS{end}.Km;                                         % Minus characteristic constant [deg]
    dataS{L}.Kp    = dataENS{L,L+1}.Kp;                                     % Plus characteristic constant [deg]
    dataS{L}.theta = 0.5*(dataS{L}.Km + dataS{L}.Kp);                       % Flow angle w.r.t horizontal [deg]
    dataS{L}.nu    = 0.5*(dataS{L}.Km - dataS{L}.Kp);                       % Prandtl-Meyer angle [deg]
    dataS{L}.M     = PM_EQUATION(dataS{L}.nu,0,g);                          % Mach number []
    dataS{L}.mu    = asind(1/dataS{L}.M);                                   % Mach angle [deg]
end

% Print the time taken to solve the straightening region
fprintf('%2.2f [sec]\n',toc);                                               % Print the time

%% FIND SHAPE OF NOZZLE [X and Y]

% Indicate that we are finding the nozzle shape, start section timer
fprintf('Solving for nozzle shape: ');                                      % Print to the command window
tic;                                                                        % Start the timer

% ==========================================
% ==== Expansion and Non-Simple Regions ====
% ==========================================
startR = 1;
for L = 2:1:numChar+1
    for R = startR:1:numChar
        if (L-1 == R)
            dataENS{R,L}.tsim = tand(0.5*((dataENS{R,L-1}.theta-dataENS{R,L-1}.mu)+...
                                       (dataENS{R,L}.theta-dataENS{R,L}.mu)));
            dataENS{R,L}.X    = dataENS{R,L-1}.X - (dataENS{R,L-1}.Y/dataENS{R,L}.tsim);
        else
            dataENS{R,L}.tsim = tand(0.5*((dataENS{R,L-1}.theta-dataENS{R,L-1}.mu)+...
                                       (dataENS{R,L}.theta-dataENS{R,L}.mu)));
            
            dataENS{R,L}.tsip = tand(0.5*((dataENS{R-1,L}.theta+dataENS{R-1,L}.mu)+...
                                       (dataENS{R,L}.theta+dataENS{R,L}.mu)));
            
            num = dataENS{R-1,L}.Y - dataENS{R,L-1}.Y + ...
                    (dataENS{R,L}.tsim*dataENS{R,L-1}.X) - ...
                        (dataENS{R,L}.tsip*dataENS{R-1,L}.X);
            den = dataENS{R,L}.tsim - dataENS{R,L}.tsip;
            
            dataENS{R,L}.X = num/den;
            dataENS{R,L}.Y = dataENS{R,L-1}.Y + ...
                            dataENS{R,L}.tsim*(dataENS{R,L}.X-dataENS{R,L-1}.X);
        end
    end
    startR = startR + 1;
end

% ==============================
% ==== Straightening Region ====
% ==============================
for L = 1:1:numChar
    if (L == 1)
        dataS{L}.tsi  = tand(0.5*(dataS{L}.theta + dataENS{numChar,1}.theta));
        dataS{L}.tsip = tand(dataENS{numChar,L+1}.theta + dataENS{numChar,L+1}.mu);
        
        num = (dataS{L}.tsi*dataENS{numChar,L}.X) - ...
              (dataS{L}.tsip*dataENS{numChar,L+1}.X) - ...
              (dataENS{numChar,L}.Y) + ...
              (dataENS{numChar,L+1}.Y);
        den = dataS{L}.tsi - dataS{L}.tsip;
        
        dataS{L}.X = num/den;
        dataS{L}.Y = dataENS{numChar,L}.Y + ...
                     (dataS{L}.tsi*dataS{L}.X) - ...
                     (dataS{L}.tsi*dataENS{numChar,L}.X);
    else
        dataS{L}.tsi  = tand(0.5*(dataS{L-1}.theta + dataS{L}.theta));
        dataS{L}.tsip = tand(dataENS{numChar,L+1}.theta + dataENS{numChar,L+1}.mu);
        
        num = (dataS{L}.tsi*dataS{L-1}.X) - ...
              (dataS{L}.tsip*dataENS{numChar,L+1}.X) - ...
              (dataS{L-1}.Y) + ...
              (dataENS{numChar,L+1}.Y);
        den = dataS{L}.tsi - dataS{L}.tsip;
        
        dataS{L}.X = num/den;
        dataS{L}.Y = dataS{L-1}.Y + dataS{L}.tsi*(dataS{L}.X - dataS{L-1}.X);
    end
end

% Print the time taken to solve the nozzle shape
fprintf('%2.2f [sec]\n',toc);                                               % Print the time

%% DISPLAY SOME RESULTS

M_exit      = dataS{numChar}.M;
MoC_A_Astar = (dataS{numChar}.Y)/(Dstar/2);
A_Astar     = A_M_RELATION(0,M_exit,g,'Sup');

fprintf('============ RESULTS ==============\n');
fprintf('Me [MoC]   : %1.5f\n',M_exit);
fprintf('Me [Q1D]   : %1.5f\n',Me_Set);
fprintf('A/A* [MoC] : %1.5f\n',MoC_A_Astar);
fprintf('A/A* [Q1D] : %1.5f\n',Ae_At_Set);

%% PLOT PATCHES OF MACH NUMBER

% Indicate that we are plotting the nozzle
fprintf('Plotting...\n');                                                   % Print to the command window

% Set up the figure
figure(1);                                                                  % Select the appropriate figure
cla; hold on; grid on;                                                      % Get ready for plotting
c = colorbar;                                                               % Turn on the colorbar object
% axis('equal');
xlim('auto');                                                               % Auto-limit the X-axis
ylim('auto');                                                               % Auto-limit the Y-axis

% Starting Constant Section
% - Where the input Mach number is M1
% - Have not reached the first characteristic
ss(1,1) = 0;                                                                % \
ss(1,2) = Dstar/2;                                                          %  |-> 
ss(1,3) = M1;                                                               % /
ss(2,1) = ss(1,1);                                                          % \
ss(2,2) = -ss(1,2);                                                         %  |-> 
ss(2,3) = ss(1,3);                                                          % /
ss(3,1) = dataENS{1,2}.X;                                                   % \
ss(3,2) = dataENS{1,2}.Y;                                                   %  |-> 
ss(3,3) = M1;                                                               % /

% Ending Constant Section
% - Where the output Mach number is Mexit
% - Have reached horizontal, fully expanded flow
se(1,1) = dataENS{numChar,numChar+1}.X;                                     % \
se(1,2) = dataENS{numChar,numChar+1}.Y;                                     %  |-> 
se(1,3) = dataENS{numChar,numChar+1}.M;                                     % /
se(2,1) = dataS{numChar}.X;                                                 % \
se(2,2) = 0;                                                                %  |-> 
se(2,3) = dataS{numChar}.M;                                                 % /
se(3,1) = dataS{numChar}.X;                                                 % \
se(3,2) = dataS{numChar}.Y;                                                 %  |-> 
se(3,3) = dataS{numChar}.M;                                                 % /

% Plot the starting and ending uniform flow patches
patch(se(:,1),se(:,2),se(:,3),'EdgeColor','none','FaceColor','interp');     % Plot the upper starting patch
patch(se(:,1),-se(:,2),se(:,3),'EdgeColor','none','FaceColor','interp');    % Plot the lower starting patch
patch(ss(:,1),ss(:,2),ss(:,3),'EdgeColor','none','FaceColor','interp');     % Plot the upper ending patch
patch(ss(:,1),-ss(:,2),ss(:,3),'EdgeColor','none','FaceColor','interp');    % Plot the lower ending patch

% Expansion and non-simple region patch plotting
sol = [];                                                                   % Initialize the solution array
for R = 1:1:numChar-1
    for L = 1:1:numChar
        if (R == L-1)                                                       % If we are on the nozzle edge
            sol(1,1) = dataENS{R,L}.X;
            sol(1,2) = dataENS{R,L}.Y;
            sol(1,3) = dataENS{R,L}.M;
            sol(2,1) = dataENS{R+1,L}.X;
            sol(2,2) = dataENS{R+1,L}.Y;
            sol(2,3) = dataENS{R+1,L}.M;
            sol(3,1) = dataENS{R+1,L+1}.X;
            sol(3,2) = dataENS{R+1,L+1}.Y;
            sol(3,3) = dataENS{R+1,L+1}.M;
        else                                                                % If we are not on the nozzle edge
            sol(1,1) = dataENS{R,L}.X;
            sol(1,2) = dataENS{R,L}.Y;
            sol(1,3) = dataENS{R,L}.M;
            sol(2,1) = dataENS{R+1,L}.X;
            sol(2,2) = dataENS{R+1,L}.Y;
            sol(2,3) = dataENS{R+1,L}.M;
            sol(3,1) = dataENS{R+1,L+1}.X;
            sol(3,2) = dataENS{R+1,L+1}.Y;
            sol(3,3) = dataENS{R+1,L+1}.M;
            sol(4,1) = dataENS{R,L+1}.X;
            sol(4,2) = dataENS{R,L+1}.Y;
            sol(4,3) = dataENS{R,L+1}.M;
        end
        patch(sol(:,1),sol(:,2),sol(:,3),'EdgeColor','none',...
                                         'FaceColor','interp');
        patch(sol(:,1),-sol(:,2),sol(:,3),'EdgeColor','none',...
                                          'FaceColor','interp');
    end
end

% Straightening region patch plotting
sol = [];
for L = 1:1:numChar
    if (L == 1)
        sol(1,1) = dataENS{numChar,L}.X;
        sol(1,2) = dataENS{numChar,L}.Y;
        sol(1,3) = dataENS{numChar,L}.M;
        sol(2,1) = dataS{L}.X;
        sol(2,2) = dataS{L}.Y;
        sol(2,3) = dataS{L}.M;
        sol(3,1) = dataENS{numChar,L+1}.X;
        sol(3,2) = dataENS{numChar,L+1}.Y;
        sol(3,3) = dataENS{numChar,L+1}.M;
    else
        sol(1,1) = dataENS{numChar,L}.X;
        sol(1,2) = dataENS{numChar,L}.Y;
        sol(1,3) = dataENS{numChar,L}.M;
        sol(2,1) = dataS{L-1}.X;
        sol(2,2) = dataS{L-1}.Y;
        sol(2,3) = dataS{L-1}.M;
        sol(3,1) = dataS{L}.X;
        sol(3,2) = dataS{L}.Y;
        sol(3,3) = dataS{L}.M;
        sol(4,1) = dataENS{numChar,L+1}.X;
        sol(4,2) = dataENS{numChar,L+1}.Y;
        sol(4,3) = dataENS{numChar,L+1}.M;
    end
    patch(sol(:,1),sol(:,2),sol(:,3),'EdgeColor','none',...
                                     'FaceColor','interp');
    patch(sol(:,1),-sol(:,2),sol(:,3),'EdgeColor','none',...
                                      'FaceColor','interp');
end

%% PLOT THE NOZZLE AND CHARACTERISTICS

% Maximum Mach number so we can set Z-axis level correctly
maxM = dataS{numChar}.M;

% ===============================
% ==== NOZZLE BOUNDING LINES ====
% ===============================

xNoz = [dataENS{1,1}.X];
yNoz = [dataENS{1,1}.Y];

% Nozzle bounding lines in expansion section
for R = 1:1:numChar-1
    xNoz = [xNoz; dataENS{R+1,1}.X];
    yNoz = [yNoz; dataENS{R+1,1}.Y];
    plot([dataENS{R,1}.X; dataENS{R+1,1}.X],...                             % Top surface
          [dataENS{R,1}.Y; dataENS{R+1,1}.Y],...
          'k-','LineWidth',3);
    plot([dataENS{R,1}.X; dataENS{R+1,1}.X],...                             % Bottom surface
          [-dataENS{R,1}.Y; -dataENS{R+1,1}.Y],...
          'k-','LineWidth',3);
end

% Nozzle bounding line between expansion and straightening sections
plot([dataENS{numChar,1}.X; dataS{1}.X],...                                 % Top surface
    [dataENS{numChar,1}.Y; dataS{1}.Y],'k-','LineWidth',3);
plot([dataENS{numChar,1}.X; dataS{1}.X],...                                 % Bottom surface
    [-dataENS{numChar,1}.Y; -dataS{1}.Y],'k-','LineWidth',3);

% Nozzle bounding lines in straightening section
xNoz = [xNoz; dataS{1}.X];
yNoz = [yNoz; dataS{1}.Y];
for L = 2:1:numChar
    xNoz = [xNoz; dataS{L}.X];
    yNoz = [yNoz; dataS{L}.Y];
    plot([dataS{L-1}.X; dataS{L}.X],...                                     % Top surface
         [dataS{L-1}.Y; dataS{L}.Y],'k-','LineWidth',3);
    plot([dataS{L-1}.X; dataS{L}.X],...                                     % Bottom surface
         [-dataS{L-1}.Y; -dataS{L}.Y],'k-','LineWidth',3);
end

% ==============================
% ==== CHARACTERISTIC LINES ====
% ==============================

% Non-simple region C- characteristics
for R = 1:1:numChar
    Xvals = [];
    Yvals = [];
    for L = 1:1:R+1
        Xvals = [Xvals; dataENS{R,L}.X];
        Yvals = [Yvals; dataENS{R,L}.Y];
    end
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

% Non-simple region C+ characteristics
for L = 2:1:numChar
    Xvals = [];
    Yvals = [];
    for R = L-1:1:numChar
        Xvals = [Xvals; dataENS{R,L}.X];
        Yvals = [Yvals; dataENS{R,L}.Y];
    end
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

% Expansion region C+ characteristics
for L = 2:1:numChar+1
    Xvals = [dataENS{numChar,L}.X; dataS{L-1}.X];
    Yvals = [dataENS{numChar,L}.Y; dataS{L-1}.Y];
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

% ==========================
% ==== EXTRANEOUS STUFF ====
% ==========================

% Nozzle centerline
plot([0; dataS{numChar}.X],[0; 0],'k-','LineWidth',3);

% Plot white dots (black outline) for intersections in non-simple region
for L = 2:1:numChar+1
    for R = 1:1:numChar
            plot(dataENS{R,L}.X,dataENS{R,L}.Y,'o',...
                  'MarkerFaceColor','w','MarkerEdgeColor','k');
    end
end

% PLOT JUST THE OUTSIDE CONTOUR
plot(xNoz,yNoz,'o','MarkerFaceColor','g','MarkerEdgeColor','k');

% Colorbar properties
caxis([M1 dataS{numChar}.M]);
ylabel(c,'Mach Number','FontSize',24);
set(c,'YTick',[M1 get(c,'YTick') dataS{numChar}.M]);

fprintf('SOLUTION FINISHED!\n');

%% LOAD CFD CSV FILE AND COMPARE

% flpth = 'C:\Users\Josh\Documents\MATLAB\Nozzle_Design\';
% flnm  = 'CSV_Data.csv';
% fid   = fopen([flpth flnm]);
% 
% dataCFD = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f',...
%                         'HeaderLines',1,'CollectOutput',1,'Delimiter',',');
% 
% % Get Mach number from symmetry line MoC
% ind = 1;
% for i = 1:1:numChar
%     for j = 1:1:numChar
%         if (i == j+1)
%             dataX(ind) = dataENS{i,j}.X;
%             dataM(ind) = dataENS{i,j}.M;
%             ind = ind + 1;
%         end
%     end
% end
% 
% figure(30);
% cla; hold on; grid on;
% plot(dataX,dataM,'k-','LineWidth',2);
% plot(dataCFD{1,1}(:,11),dataCFD{1,1}(:,8),'ro');

%% OUTPUT GEO MESH FILE FOR GMSH

if (outputMesh == 1)

% Need to write the following
% - Point array
% - Line array
% - Line Loop
% - Plane Surface
% - Physical Line (Inlet, Outlet, Symmetry, Nozzle)

% Nozzle boundary points
% - Get rid of repeated values (from point-corner for example)
xNoz = unique(xNoz);                                                        % Nozzle boundary X-points
yNoz = unique(yNoz);                                                        % Nozzle boundary Y-points

% Number of nozzle points and lines
numNozPts = length(xNoz);                                                   % Number of nozzle points
numNozLns = length(xNoz)-1;                                                 % Number of nozzle lines

% Clear any old variables
clearvars Point;

% ===== FILE PROPERTIES =====
fid = fopen('Output_Mesh.geo','w');

% ===== POINTS =====

% Nozzle bounday
Point(:,1) = xNoz;
Point(:,2) = yNoz;

% Additional points
addPts = [xNoz(end) 0;
          0         0;
          xNoz(1)   yNoz(1)];
Point  = [Point; addPts];

Point(:,3) = 0.0;
Point(:,4) = 1.0;

% Number of points
numPts = size(Point,1);

% Write points to a file
for i = 1:1:numPts
    fprintf(fid,'Point(%i) = {%g, %g, %g, %g};\r\n',...
                    i,Point(i,1),Point(i,2),Point(i,3),Point(i,4));
end

% ===== LINES =====

% Number of lines
numLns = numPts;

% Write lines to file
for i = 1:1:numPts
    if (i ~= numLns)
        fprintf(fid,'Line(%i) = {%i, %i};\r\n',i,i,i+1);
    else
        fprintf(fid,'Line(%i) = {%i, %i};\r\n',i,i,1);
    end
end

% ===== LINE LOOP =====

% Create string with all lines
lineStr = '1';
for i = 2:1:numLns
    lineStr = [lineStr ', ' num2str(i)];
end

% Write to file
fprintf(fid,'Line Loop(%i) = {%s};\r\n',numPts+1,lineStr);

% ===== PLANE SURFACE =====

fprintf(fid,'Plane Surface(%i) = {%i};\r\n',numPts+2,numPts+1);

% ===== PHYSICAL LINE =====

phyLineNozzle = '1';
for i = 2:1:numNozLns
    phyLineNozzle = [phyLineNozzle ', ' num2str(i)];
end

fprintf(fid,'Physical Line("Outlet") = {%i};\r\n',numLns-3);
fprintf(fid,'Physical Line("Symmetry") = {%i};\r\n',numLns-2);
fprintf(fid,'Physical Line("Inlet") = {%i};\r\n',numLns-1);
fprintf(fid,'Physical Line("Nozzle") = {%s};\r\n',phyLineNozzle);

% ===== CHARACTERISTIC LENGTH =====

charVal    = 0.2;
charLength = '1';
for i = 2:1:numPts
    charLength = [charLength ', ' num2str(i)];
end

fprintf(fid,'Characteristic Length {%s} = %f;\n',charLength,charVal);

end

fclose('all');






