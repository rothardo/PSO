clc, clear;


costFcnLower =@(x) CostFunction(x);
nvars = 2;      
LB = [-1 -1];   
UB = [ 2  2];  

opt = optimoptions( ...
    'particleswarm','SwarmSize',50,'MaxIterations',200*nvars, ...
    'MaxStallIterations',20, 'HybridFcn',@fmincon, ...
    'UseParallel',false,'Display','iter','PlotFcn', ...
    @pswplotbestf); 
[xOptimal,zOptimal] = particleswarm(costFcnLower,nvars,LB,UB,opt);
function Z = CostFunction(x)

z1 = sin(x(1));

z2 = sec(x(2));

Z = z1+z2;
end
