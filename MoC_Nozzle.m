% Method of Characteristics - Nozzle Design
% Written by: JoshTheEngineer
% Started: 06/24/16
% Updated: 06/24/16 - Started code
%                   - Works as expected, but no plotting yet
%          06/25/16 - Plotting works now
%          06/27/16 - Added comments

clear;
clc;

%% INPUTS

% Shoulder angles
theta_max = 20;                                                             % Maximum throat expansion angle [deg]
Dstar     = 2;                                                              % Diameter [length]
geom      = 2;                                                              % [1] Point, [2] Circle
radius    = 1;                                                              % Radius of throat [length]
numChar   = 15;                                                             % Number of characteristic lines to use
theta     = linspace(0,theta_max,numChar)';                                 % Theta array for each characteristic

% Throat Mach number
M1 = 1.05;

% Specific heat ratio
g = 1.4;

%% INITIAL SETUP OF KNOWNS

% Initialize solution variables
dataENS  = cell(numChar,numChar+1);                                            % Expansion and non-simple regions
dataS = cell(numChar,1);                                                    % Straightening region

% ====================
% ==== EXPANSION =====
% ====================
for R = 1:1:numChar                                                         % Loop through all (-) characteristics
    dataENS{R,1}.M = 0;                                                     % Set all Mach numbers to zero
    if (R == 1)                                                             % For the first point
        dataENS{R,1}.M = M1;                                                    % Set the Mach number from the throat M
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
        dataENS{R,1}.X = 0;                                                     % All X-values are at the throat
        dataENS{R,1}.Y = Dstar/2;                                               % All Y-values are at half the throat diameter
    elseif (geom == 2)                                                      % Finite-radius corner
        dataENS{1,1}.X = 0;                                                     % First X-value is at the throat
        dataENS{1,1}.Y = Dstar/2;                                               % First Y-value is at half the throat diameter
        
        for i = 2:1:numChar                                                     % For the rest of the starting line characteristics
            dx = radius*sind(theta(i));                                             % Change in X depends on radius and angle
            dy = radius - radius*cosd(theta(i));                                    % Change in Y depends on radius and angle
            
            dataENS{i,1}.X = dataENS{1,1}.X + dx;                                   % Add the X change to the first X point
            dataENS{i,1}.Y = dataENS{1,1}.Y + dy;                                   % Add the Y change to the first Y point
        end
    end
end

% =====================
% ==== NON-SIMPLE =====
% =====================
for R = 1:1:numChar                                                         % Loop over all (-) characteristics
    for L = 2:1:numChar+1                                                   % Loop over all (+) characteristics
        dataENS{R,L}.M     = 0;                                                 % Mach number
        dataENS{R,L}.theta = 0;                                                 % Flow angle w.r.t horizontal [deg]
        dataENS{R,L}.nu    = 0;                                                 % Prandtl-Meyer angle [deg]
        dataENS{R,L}.mu    = 0;                                                 % Mach angle [deg]
        dataENS{R,L}.Kp    = 0;                                                 % Plus characteristic constant [deg]
        dataENS{R,L}.Km    = 0;                                                 % Minus characteristic constant [deg]
        dataENS{R,L}.X     = 0;                                                 % Point X location
        dataENS{R,L}.Y     = 0;                                                 % Point Y location
        dataENS{R,L}.dydx  = 0;                                                 % Point slope
        dataENS{R,L}.tsip  = 0;                                                 % Convenient avg parameter for (+) characteristic
        dataENS{R,L}.tsim  = 0;                                                 % Convenient avg parameter for (-) characteristic
    end
end

% ========================
% ==== STRAIGHTENING =====
% ========================
for L = 1:1:numChar                                                         % Loop through all (+) characteristics
    dataS{L}.M     = 0;                                                         % Mach number []
    dataS{L}.theta = 0;                                                         % Flow angle w.r.t. horizontal[deg]
    dataS{L}.nu    = 0;                                                         % Prandtl-Meyer angle [deg]
    dataS{L}.mu    = 0;                                                         % Mach angle [deg]
    dataS{L}.Kp    = 0;                                                         % Plus characteristic constant [deg]
    dataS{L}.Km    = 0;                                                         % Minus characteristic constant [deg]
    dataS{L}.X     = 0;                                                         % Point X location
    dataS{L}.Y     = 0;                                                         % Point Y location
    dataS{L}.dydx  = 0;                                                         % Point slope
    dataS{L}.tsi   = 0;                                                         % Convenient nozzle contour parameter
    dataS{L}.tsip  = 0;                                                         % Convenient avg parameter for (+) characteristic
    dataS{L}.tsim  = 0;                                                         % Convenient avg parameter for (-) characteristic
end

%% EXPANSION REGION

% Indicate that we are solving expansion region, start section timer
fprintf('Solving expansion region: ');                                      % Print to the command window
tic;                                                                        % Start the timer

% Starting point Prandtl-Meyer angle
dataENS{1,1}.nu = PM_EQUATION(0,dataENS{1,1}.M,g);                          % PM angle from the throat Mach number [deg]

% All throat K+ characteristic constants
for R = 1:1:numChar                                                         % Loop over all the (-) characteristics
    dataENS{R,1}.Kp = dataENS{1,1}.theta - dataENS{1,1}.nu;                     % All the same - coming from throat [deg]
end

% PM angles, Mach numbers, Mach angles, and K- characterstic constants
for R = 1:1:numChar                                                         % Loop over all the (-) characteristics
    if (R ~= 1)                                                             % If we are not on the first characteristic (values already defined)
        dataENS{R,1}.nu = dataENS{R,1}.theta - dataENS{R,1}.Kp;                 % Prandtl-Meyer angle (don't need first characteristic)
        dataENS{R,1}.M  = PM_EQUATION(dataENS{R,1}.nu,0,g);                  	% Mach number (don't need first characteristic)
    end
    dataENS{R,1}.mu = asind(1/dataENS{R,1}.M);                                	% Mach angle [deg]
    dataENS{R,1}.Km = dataENS{R,1}.theta + dataENS{R,1}.nu;                   	% Plus characteristic constant [deg]
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
            dataENS{R,L}.nu = dataENS{R,1}.Km;                              	% Prandtl-Meyer angle from (-) constant [deg]
            dataENS{R,L}.M  = PM_EQUATION(dataENS{R,L}.nu,0,g);               	% Mach number from PM equation useing nu as input []
            dataENS{R,L}.mu = asind(1/dataENS{R,L}.M);                        	% Mach angle [deg]
            dataENS{R,L}.Kp = dataENS{R,L}.theta - dataENS{R,L}.nu;         	% Plus characteristic constant [deg]
            dataENS{R,L}.Km = dataENS{R,L}.theta + dataENS{R,L}.nu;            	% Minus characteristic constant [deg]
        else                                                                % For all other (+) characteristics
            dataENS{R,L}.Kp    = dataENS{R-1,L}.Kp;                            	% Plus characteristic constant [deg]
            dataENS{R,L}.Km    = dataENS{R,L-1}.Km;                           	% Minus characteristic constant [deg]
            dataENS{R,L}.theta = 0.5*(dataENS{R,L}.Km + dataENS{R,L}.Kp);    	% Angle w.r.t. horizontal [deg]
            dataENS{R,L}.nu    = 0.5*(dataENS{R,L}.Km - dataENS{R,L}.Kp);     	% Prandtl-Meyer angle [deg]
            dataENS{R,L}.M     = PM_EQUATION(dataENS{R,L}.nu,0,g);             	% Mach number []
            dataENS{R,L}.mu    = asind(1/dataENS{R,L}.M);                      	% Mach angle [deg]
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
            dataENS{R,L}.X = dataENS{R,L-1}.X - (dataENS{R,L-1}.Y/dataENS{R,L}.tsim);
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
A_Astar     = sqrt((1/M_exit^2)*...
               (((2/(g+1)*(1+((g-1)/2)*M_exit^2))^((g+1)/(g-1)))));

fprintf('MoC Mexit : %1.2f\n\n',M_exit);
fprintf('MoC A/A*  : %1.4f\n',MoC_A_Astar);
fprintf('Thry A/A* : %1.5f\n',A_Astar);

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
ss(3,1) = dataENS{1,2}.X;                                                      % \
ss(3,2) = dataENS{1,2}.Y;                                                      %  |-> 
ss(3,3) = M1;                                                               % /

% Ending Constant Section
% - Where the output Mach number is Mexit
% - Have reached horizontal, fully expanded flow
se(1,1) = dataENS{numChar,numChar+1}.X;                                        % \
se(1,2) = dataENS{numChar,numChar+1}.Y;                                        %  |-> 
se(1,3) = dataENS{numChar,numChar+1}.M;                                        % /
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

% Nozzle bounding lines in expansion section
for R = 1:1:numChar-1
    plot([dataENS{R,1}.X; dataENS{R+1,1}.X],...
          [dataENS{R,1}.Y; dataENS{R+1,1}.Y],...
          'k-','LineWidth',3);
    plot([dataENS{R,1}.X; dataENS{R+1,1}.X],...
          [-dataENS{R,1}.Y; -dataENS{R+1,1}.Y],...
          'k-','LineWidth',3);
end

% Nozzle bounding line between expansion and straightening sections
plot([dataENS{numChar,1}.X; dataS{1}.X],...
    [dataENS{numChar,1}.Y; dataS{1}.Y],'k-','LineWidth',3);
plot([dataENS{numChar,1}.X; dataS{1}.X],...
    [-dataENS{numChar,1}.Y; -dataS{1}.Y],'k-','LineWidth',3);

% Nozzle bounding lines in straightening section
for L = 2:1:numChar
    plot([dataS{L-1}.X; dataS{L}.X],...
         [dataS{L-1}.Y; dataS{L}.Y],'k-','LineWidth',3);
    plot([dataS{L-1}.X; dataS{L}.X],...
         [-dataS{L-1}.Y; -dataS{L}.Y],'k-','LineWidth',3);
end

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

% Nozzle centerline
plot([0; dataS{numChar}.X],[0; 0],'k-','LineWidth',3);

% Plot white dots (black outline) for intersections in non-simple region
for L = 2:1:numChar+1
    for R = 1:1:numChar
            plot(dataENS{R,L}.X,dataENS{R,L}.Y,'o',...
                  'MarkerFaceColor','w','MarkerEdgeColor','k');
    end
end

% Colorbar properties
caxis([M1 dataS{numChar}.M]);
ylabel(c,'Mach Number','FontSize',24);
set(c,'YTick',[M1 get(c,'YTick') dataS{numChar}.M]);

fprintf('SOLUTION FINISHED!');

