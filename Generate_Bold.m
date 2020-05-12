function [b] = Generate_Bold(T,u,dt)

% The Hemodynamic model with one simplified neural activity
% T       : total time (s)
% u       : synaptic activity
% dt      : Step size for Euler's Method

t0  = (0:dt:T)';
n_t = length(t0);
t_min = 20;
n_min = round(t_min/dt);

% BOLD model parameters

% eps    = 1.5*0.54;    % efficacy
% taus   = 1.54;    % time unit (s)
% tauf   = 2.46;    % time unit (s)
% tauo   = 2*0.98;    % mean transit time (s)
% alpha  = 0.33;    % a stiffness exponent
% Eo     = 0.34;    % resting oxygen extraction fraction

eps    = 1.5*0.54;    % efficacy
taus   = 0.65;    % time unit (s)
tauf   = 0.41;    % time unit (s)
tauo   = 0.98;    % mean transit time (s)
alpha  = 0.32;    % a stiffness exponent
Eo     = 0.34;    % resting oxygen extraction fraction



itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;

vo     = 0.02;
k1     = 7*Eo; 
k2     = 2; 
k3     = 2*Eo-0.2;

% Initial conditions

x0  = [0 1 1 1];


% Euler method

x      = zeros(n_t,4);
x(1,:) = x0;
for n = 1:n_t-1
    x(n+1,1) = x(n,1) + dt*( eps*u(n)-itaus*x(n,1)-itauf*(x(n,2)-1) );
    x(n+1,2) = x(n,2) + dt*x(n,1);
    x(n+1,3) = x(n,3) + dt*itauo*(x(n,2)-x(n,3)^ialpha);
    x(n+1,4) = x(n,4) + dt*itauo*(x(n,2)*(1-(1-Eo)^(1/x(n,2)))/Eo - (x(n,3)^ialpha)*x(n,4)/x(n,3));
    %[x(n+1,1) x(n+1,2) x(n+1,3) x(n+1,4)]
end

v  = x(n_min:end,3);
q  = x(n_min:end,4);
b  = 100/Eo*vo*( k1*(1-q) + k2*(1-q./v) + k3*(1-v) );