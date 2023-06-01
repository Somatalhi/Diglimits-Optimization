% (MS Thesis)
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dig Limits Optimization Using Binary Integer Linear
% Programming Method In Open Pit Mines 
% Author: Hussam N. Altalhi
% Advisor: Awuah-Offei, Kwame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------
%  Mathematical Model:(BILP)
%
%  Maximize:Profit 
%
%  Subjected to:4 constraints ((3x) Shape constraints + (1x) Destination
%   Constraint) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sol,fval,exitflag, t_build, t_solve] = diglimitsv6(d, alpha, beta)
%--------------------------------------------------
% % The name of the function is 'diglimitsv6'
% This function has 3 inputs and 5 outputs.
%
% Inputs:
% d (Number of Destinations)
% alpha (Minimum Number of Blocks (along i-direction)
% beta (Minimuum Number of Blocks (along j-direction)
%
% Outputs: 
% sol (The values of the decision variables)
% fval (The value of the objective function)
% exitflag (The optimal solution found)
% t_build (The time to formulate the problem)
% t_solve (The time to find the optimal solution)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************************
% requires: The user has to make sure blocks at the boundary are waste
%           (First row & First column are waste blocks).
%
%************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the data for the optimization problem
%--------------------------------------------------
% Block Model and Input Parameter file
finput = -1;
while finput < 0        % repeat until valid file name entered
    [inputfile, inputpath] = uigetfile('*.xlsx',...
        'Select input file');  % get GUI file selector
    finput = fopen([inputpath inputfile],'r');    % read only
    if finput < 0, disp('! File not found or file sharing problem !'), end
end

opts                = detectImportOptions([inputpath inputfile]);
opts.MissingRule    = "omitrow";
data                = readtable([inputpath inputfile],opts);

fclose(finput); % Close file

I = data.I;
J = data.J;
n = max(I); %number of blocks along i-direction
w = max(J); %number of blocks along j-direction

% Rewrite the the data into tables 
v = zeros(n,w,d);
[a,~] = size(data);
 
for i = 1:a
    v(data.I(i),data.J(i),:) = [data.V1(i) data.V2(i) data.V3(i) data.V4(i) data.V5(i) data.V6(i)];
end

tstart = tic; % initialize time to formulate problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the Optimization Problem
%--------------------------------------------------
prob  = optimproblem('ObjectiveSense','maximize');
y       = optimvar('y',n,w,d,'Type','integer','LowerBound',0,'UpperBound',1);
z1      = optimvar('z1',n,w,d,'Type','integer','LowerBound',0,'UpperBound',1);
z2      = optimvar('z2',n,w,d,'Type','integer','LowerBound',0,'UpperBound',1);

% Objective function: Maximize the profit.
%--------------------------------------------------
prob.Objective    = 0;
for i = 1:d
    prob.Objective = prob.Objective + sum(sum(v(:,:,i).* y(:,:,i))); 
end


% Constraint1: Each block has to be sent at most to one destination.
%--------------------------------------------------
prob.Constraints.cons1 = sum(y,3) <= 1;

% Constraint2: Block(ij) must be mined to be the left most point.
%--------------------------------------------------
prob.Constraints.cons2v=optimconstr(n,w,d);
prob.Constraints.cons2h=optimconstr(n,w,d);
prob.Constraints.cons2v = z1 - y <= 0;
prob.Constraints.cons2h = z2 - y <= 0;

% Constraint3: Assuring that the added waste blocks at the boundary are not
% being mined. 
%--------------------------------------------------

prob.Constraints.cons4v=optimconstr(n,w,d);
prob.Constraints.cons4v = y(1,:,:) <= 0;
prob.Constraints.cons4h = y(:,1,:) <= 0;

  
% Constraint4: if block (ij) is mined but (ij-1) is not then block(ij) is the
% leftmost block (Contiguity) && Constraint5: Mininmum Width of Mined blocks (Equipment
% ability of digging). (shape constraints)  
%--------------------------------------------------

%Vertical Constraints.
prob.Constraints.cons5v=optimconstr(n,w,d);
prob.Constraints.cons6v=optimconstr(n,w,d);

for i=1:n
    for j=1:w
        for k = 1:d
            if i > 1; prob.Constraints.cons5v(i,j,k) = y(i,j,k) - y(i-1,j,k) <= z1(i,j,k); end
            prob.Constraints.cons6v(i,j,k) = sum(z1(max(i-alpha,1):i,j,k)) <= y(i,j,k);
        end
    end
end

%Horizontal Constraints.
prob.Constraints.cons5h=optimconstr(n,w,d);
prob.Constraints.cons6h=optimconstr(n,w,d);

for i=1:n
    for j=1:w
        for k = 1:d
            if j > 1; prob.Constraints.cons5h(i,j,k) = y(i,j,k) - y(i,j-1,k) <= z2(i,j,k); end
            prob.Constraints.cons6h(i,j,k) = sum(z2(i,max(j-beta,1):j,k)) <= y(i,j,k);
        end
    end
end


t_build = toc(tstart);

tstart = tic;

options                 = optimoptions('intlinprog');
[sol, fval,exitflag]    = solve(prob,'options',options);

t_solve = toc(tstart);

% Removing Waste blocks at the boundary
sol.y(:,1,:)=[];
sol.y(1,:,:)=[];
sol.z1(:,1,:)=[];
sol.z1(1,:,:)=[];
sol.z2(:,1,:)=[];
sol.z2(1,:,:)=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of Optimization Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
