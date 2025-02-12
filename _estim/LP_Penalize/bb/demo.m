% Smooth Local Projections Demo
% R. Barnichon and C. Brownlees, 04/2018

clc
clear 

%% Load data
data = csvread('data.csv',1,1);

T = size(data,1);
P  = 4; % number of lags used in LP for controls

%% Estimating the IR of GDP to a monetary shock

H_min = 1; % start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_max = 20; 

y  = data(:,1); % endogenous variable
x  = data(:,3); % endoegnous variable related to the shock 
w  = [ data(:,1:2) , lagmatrix( data , 1:P ) ]; % control variables (contemporaneous vars, lagged vars)
newData = cat(2, y, x, w)

% Remove missings from data
newData(any(isnan(newData), 2), :) = [];

% Re-declare variables after removing missings
y  = newData(:,1); % endogenous variable
x  = newData(:,2); % endoegnous variable related to the shock
w = newData(:,3:size(newData,2)); % control variables and lags

lp = locproj(y,x,w,H_min,H_max,'reg'); % IR from (standard) Local Projection

r = 2; %(r-1)=order of the limit polynomial (so r=2 implies the IR is shrunk towards a line )
lambda = 100; % value of penalty

slp    = locproj(y,x,w,H_min,H_max,'smooth',r,lambda); %IR from Smooth Local Projection
slp_lim= locproj(y,x,w,H_min,H_max,'smooth',r,1e10); % Limit IR in Smooth Local Projection

figure(1)
hold on,
plot( 0:H_max , [ lp.IR slp.IR slp_lim.IR] )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])
legend('IR_{lp}','IR_{slp}','IR_{slp,max pen}','Location','Best')

%% Cross-Validation Choice of Lambda

slp = locproj(y,x,w,H_min,H_max,'smooth',r,0.01);

lambda = [ 1:0.5:10] * T;
slp    = locproj_cv(slp,5,lambda);

figure(2)
plot( lambda , slp.rss , '-o' )

lambda_opt = lambda( min( slp.rss ) == slp.rss );

%% Confidence Intervals

lp = locproj(y,x,w,H_min,H_max,'reg'); % IR from (regular) Local Projection
lp = locproj_conf(lp,H_max); % it takes a bit too run! please be patient

figure(3)
hold on,
plot( 0:H_max , lp.IR   , 'r' , 'LineWidth' , 2 )
plot( 0:H_max , lp.conf , 'r' )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])

r      = 2;
slp    = locproj(y,x,w,H_min,H_max,'smooth',r,lambda_opt); 
slp    = locproj_conf(slp,H_max,lambda_opt/2);

figure(4)
hold on,
plot( 0:H_max , slp.IR   , 'r' , 'LineWidth' , 2 )
plot( 0:H_max , slp.conf , 'r' )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])
