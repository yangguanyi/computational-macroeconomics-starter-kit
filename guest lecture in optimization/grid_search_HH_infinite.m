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
    
    % same as the two period model algorithm, our purpose is to find the
    % choice of K' that makes the life time value the highest, and we store
    % the value of K' in decision rule G. 
    
    % Difference is that instead of computing directly the future (second
    % period) value, we assume it a starting value 0, and taken it as given
    % to find K' to maximize Tv = u(c) + beta*v with the assumed v. 
    
    % At the end of each iteration, we update our assumption of V (future
    % value), since we know in steady state (stationary equilibrium), Tv =
    % V. 
    
    % We successfully found our solution when Tv = V. This is done through
    % updating the "distance" our stopping criteria. 
    
    Tv0 = zeros(snum, snum); % Tv0 records all the possible candidates of lifetime value currently, given all future possible choices of s'
    c0 = Tv0;
    for i = 1:snum  % current s
        for j = 1:snum  % future s'
            c0(i,j) =  w + (1+r)*sgrid(i) - sgrid(j);
            if c0(i,j) <0
                c0(i,j)=1e-9;  % c cannot be negative
            end
            Tv0(i,j) = log(c0(i,j))+ beta*V(j);            
        end
    end
    
    % choose the s' that maximizes lifetime value
    for i = 1:snum
        [Tv(i),loc] = max(Tv0(i,:));
        s_rule(i) = sgrid(loc);
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
title ( ' the value function (grid search)' )

subplot(212)
plot(sgrid, sgrid)
hold on
plot(s_rule, sgrid)
title ( ' the decision rule ' )

saveas(gcf,'grid search.png')