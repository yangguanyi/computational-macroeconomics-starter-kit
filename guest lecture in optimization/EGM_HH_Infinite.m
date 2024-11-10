clear; close all; clc;
tic
% parameter values:
beta = 0.96;
r = 0.01; 
w = 1;


% declare the state space for savings s (exogenous grid):
sl = 5;
sh = 10;
snum = 500;
sgrid = linspace(sl, sh, snum);


% Initiate current period value: Tv, future period value: V;
Tv = zeros(size(sgrid));
V = Tv; 
% Initiate savings decision rule: 
s_rule = Tv;

% Method of successive approximation to find Tv = V: 

% define stopping rule: 
precision = 1e-5;
distance = 2*precision;
iteration = 0; % set counter to see how many iterations do we need to converge

while distance > precision
    

Dev = zeros(size(sgrid));
segrid = Dev;

for i = 1:snum 
    if i == snum
	Dev(i) = (V(i) - V(i-1))/(sgrid(i) - sgrid(i-1)); 
    else 
	Dev(i) = (V(i+1) - V(i))/(sgrid(i+1) - sgrid(i)); 
    end 

    dV_f = beta*Dev(i);
    c_temp = 1/dV_f;    %%%%% Euler equation 
    segrid(i) = c_temp + sgrid(i) - w;   %%%%% budget constraint to find current savings, which is our endogenous grid.  
end 

%%%%%% now we interpolate the exogenous grid on the endogenous grid to find the savings rule: 
for i = 1:snum 
    sval = sgrid(i);  %%%% for each s' (based on exogenous grid sgrid)

    if sval < segrid(1) 
	s_rule_val = sgrid(1);  %%%%% imposing borrowing limit
    else
	%%%%% find the position of sval on the endogenous grid
	unit = ones(size(segrid));
	resource = (sval+0.00000001)*unit - segrid;
	sindex = sum(resource > 0);
	if sindex <= 0
   	   sindex = 1;
	end

        weight = segrid(sindex+1) - sval;
        weight = weight/(segrid(sindex+1) - segrid(sindex));
        weight = min(weight, 1.0);
        weight = max(weight, 0.0);

        s_rule_val = weight*sgrid(sindex) + (1.0 - weight)*sgrid(sindex+1);

        if s_rule_val < sgrid(1)
           s_rule_val = sgrid(1);
        end
    end
    s_rule(i) = s_rule_val;
end 


%%%%% so far, we have used the EGM to find the savings decision rule s_rule. Next, we solve/update the future value based on the savings rule s_rule and the exogenous grid sgrid. 


for i = 1:snum

    sval = sgrid(i);  %%%% now this is current savings, rather than future savings 

    %%%%%% if you like, you can calculate all current period flow values
    yval = (1+r)*sval + w; %%% cash on hand
    cval = max(yval - s_rule(i), 1e-10); 
    
    u = log(cval);
	
    %%%% importantly, interpolate future value associated with s_rule and the exogenous grid sgrid: 
    
    unit = ones(size(sgrid));
    resource = (s_rule(i)+0.00000001)*unit - sgrid;
    sindex = sum(resource > 0);
    if sindex <= 0
      sindex = 1;
    end    
    
    if(sindex < snum)
      weight = sgrid(sindex+1) - s_rule(i);
      weight = weight/(sgrid(sindex+1) - sgrid(sindex));
    else
      weight = 0.0;
      sindex = sindex - 1;
    end    
   
    valf = V(sindex)*weight + V(sindex+1)*(1.0 - weight);
    Tv(i) = u + beta*valf;
end 

    
    % successive approximation -- update the future value V, using what we
    % have calculated Tv. 
    distance = max(max((abs(Tv - V))));   % We know from theory that Tv = v. "distance" is to update how good our initial guess of "v" is, given our optimization. 
    V = Tv;  % we update future "v" using our optimization of "Tv". When this current iteration is done, it will restart from the first line of "while" loop, 
             % with our initial guess of "v" coming from here (the optimal
             % value from previous iteration). We then take "v" as given,
             % and optimize again, and update again; until "Tv = v". 
    iteration = iteration + 1;   % update the counter. 
    
    % this is how I want to display the result. Telling me how many times
    % it has been iterated upon; and how close our guess is each time. 
    s = sprintf ( ' iteration %4d    ||Tv-v|| = %8.6f ', iteration, distance);
    disp(s)    
        
end

toc

figure

subplot(211)
plot(sgrid, V)
hold on
title ( ' the value function (EGM)' )

subplot(212)
plot(sgrid, sgrid)
hold on
plot(s_rule, sgrid)
title ( ' the decision rule ' )

saveas(gcf,'egm.png')