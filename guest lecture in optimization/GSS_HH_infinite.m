clear; close all; clc;
tic
% parameter values:
beta = 0.96;
r = 0.01;
w = 1;


% declare the state space for savings s:
sl = 5;
sh = 10;
snum = 500;
sgrid = linspace(sl, sh, snum);


% Initiate current period value: Tv, future period value: V;
Tv = zeros(size(sgrid));
V = Tv; 
% Initiate savings decision rule: 
s_rule = Tv;

% Method of successive approximation to find Tv = v: 

% define stopping rule: 
precision = 1e-5;
distance = 2*precision;
iteration = 0; % set counter to see how many iterations do we need to converge

while distance > precision
    
    for i = 1:snum
	    sval = sgrid(i);
        [s_rule(i), Tv(i)] = f_gss(w,sval,sgrid,r,beta,V,snum);
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
title ( ' the value function (GSS)' )

subplot(212)
plot(sgrid, sgrid)
hold on
plot(s_rule, sgrid)
title ( ' the decision rule ' )

saveas(gcf,'GSS.png')