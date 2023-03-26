%% ThreeBody_Main.m
% Author: Rishi Mayekar
% Course: ASTR 2030 "Black Holes"
% Purpose: To walk throught the process of simulating a three-body problem
%          as commonly viewed to occur in space.

clc; clear all; close all;

%% Constants and Initializations

% Constants
const.Msun = 2e30; % kg
const.G = 1; 
const.M = [1,1,1]; 
const.AU = 1.5e11; 

% initialize system
year2Sec = 365*24*3600;
const.N = 3;
x = [-0.97000436,0.,0.97000436]; %x1 = -x3, x2 = 0
y = [0.24208753,0.,-0.24208753]; %y1 = -y3, y2 = 0
vx = [0.4662036850,-0.933240737,0.4662036850]; %v1x = v3x
vy = [0.4323657300,-0.86473146,0.4323657300]; %v1y = v3y
X0 = [x, y, vx, vy];

year2Sec = 365*24*3600;
tspan = linspace(0, 6.28, 1001);


%% Function call and data handling

[T_out, X_out] = ode45(@(T_out, X_out) ThreebFun(X_out, const), tspan, X0);

figure(1);
hold on; grid on;
xlim([min(min(X_out(:,1:3))), max(max(X_out(:,1:3)))]);
ylim([min(min(X_out(:,4:6))), max(max(X_out(:,4:6)))]);

colors = ['r','k','b'];

for j = 1:length(tspan)
    for i = 1:const.N
        scatter(X_out(j,i), X_out(j,i+3), '.', colors(i))
    end
    pause(0.02)
end



%% 3 body funcition

function sdot = ThreebFun(Xstate, const)
    
    % Pulling useful constants
    G = const.G;
    M = const.M;
    N = const.N;

    x = Xstate(1:N); y = Xstate(N+1:2*N);
    u = Xstate(2*N+1:3*N); v = Xstate(3*N+1:4*N);
    
    % Change in pos is vel
    xdot = u;
    ydot = v;

    % Change in vel is acceleration
    udot = zeros(N,1); vdot = zeros(N,1);

    for i = (1:N)
        for j = (1:N)
            if i ~= j
                r_ij = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);

                udot_j = -G*M(j)*(x(i)-x(j))/r_ij^3;
                vdot_j = -G*M(j)*(y(i)-y(j))/r_ij^3;

                udot(i) = udot(i) + udot_j;
                vdot(i) = vdot(i) + vdot_j;
            end
        end
        
    end
    
    sdot = [xdot; ydot; udot; vdot];
end



