clc; clear; close all;
h= 0.1; t(1)=0; tfinal=300; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h); % Inputs
S(1) = 2.67e7;
V(1) = 9.61e6;
E(1) = 7;
I(1) = 0;
H(1) = 0;
R(1) = 0;
% Total Population
total_pop = S(1) + V(1) + E(1) + I(1) + H(1) + R(1);
alpha=0.70; % this is fractional order
tau=1; % this is fractal dimension
sigma1 = 0.114;   % Recruitment rate of individuals into the population
sigma2 = 0.000256; % Natural death rate
sigma3 = 0.958;   % Average effective contact rate
sigma4 = 0.114;   % Vaccination rate
sigma5 = 1.00;    % Vaccine inefficacy
sigma6 = 0.0015;  % Average hospitalization rate
sigma7 = 1.76;     % Recovery rate for infected individuals
sigma8 = 0.34;    % Recovery rate for exposed individuals
sigma9 = 0.092;    % Recovery rate for hospitalized individuals
sigma10 = 1 / 1.6;% Average latent or incubation rate (days^-1)
sigma11 = 1 / 50; % Rate at which individuals lose immunity (days^-1)
% The Model (3.4) in the Research Paper
f1 = @(t,S,V,E,I,H,R) sigma1 + sigma11 *(V + R) -sigma3 *(E+I).*S / total_pop - (sigma4+sigma2)*S;
f2 = @(t,S,V,E,I,H,R) sigma4*S-sigma5 * sigma3 * (E+I).*V / total_pop-(sigma11 + sigma2)*V;
f3 = @(t,S,V,E,I,H,R) sigma3 * (E+I).*S / total_pop - (sigma10 + sigma2)*E;
f4 = @(t,S,V,E,I,H,R) sigma10*E - (sigma8 + sigma2)*I;
f5 = @(t,S,V,E,I,H,R) sigma8 * I - (sigma6 + sigma9 + sigma2)*H;
f6 = @(t,S,V,E,I,H,R) sigma6*H-(sigma7 + sigma2)*R;
% Algorithm of the Caputo Fractal-Fractional starts
for n=1:N
    j=2:n;
    S(n+1)=S(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f1(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f1(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    
    V(n+1)=V(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f1(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f1(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    
    E(n+1)=E(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f1(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f1(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    
    I(n+1)=I(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f2(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f2(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    
    H(n+1)=H(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f2(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f2(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    
    R(n+1)=R(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
        f2(t(j),S(j),V(j),E(j),I(j),H(j),R(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f2(t(j-1),S(j-1),V(j-1),E(j-1),I(j-1),H(j-1),R(j-1)).*t(j-1).^(tau-1));
    t(n+1)=t(n)+h;
end
figure (1)
plot(t, S, 'r-', 'LineWidth', 3);