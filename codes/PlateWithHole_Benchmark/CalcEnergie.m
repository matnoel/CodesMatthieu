% Script qui calcul l'energie au bord du perçage pour les différents splits

clc, close all, clear all

%% Data

verif=false;

E = 210; %Pa 210e9"
v = 0.3; % ]-1;0.5[

mu = E/(2*(1+v));
lambda = v*E/((1+v)*(1-2*v));
K=E/(3*(1-2*v)); %(3*lambda+2*mu)/3

Sig = 10000; %Pa

Sig_1 = zeros(3,3);
Sig_1(2,2) = -3*Sig;

Sig_2 = zeros(3,3);
Sig_2(1,1) = Sig;

I=eye(3);

Eps_1 = (1+v)/E*Sig_1 - v/E*trace(Sig_1)*I;
Eps_2 = (1+v)/E*Sig_2 - v/E*trace(Sig_2)*I;

%% Operators

%det
detEps_1=det(Eps_1);
detEps_2=det(Eps_2);

%trace
trEps_1 = trace(Eps_1); 
trEps_2 = trace(Eps_2);

%trace p m
trEpsP_1 = (trEps_1+abs(trEps_1))/2;
trEpsP_2 = (trEps_2+abs(trEps_2))/2;

trEpsM_1 = (trEps_1-abs(trEps_1))/2;
trEpsM_2 = (trEps_2-abs(trEps_2))/2;

%deviateur
EpsDev_1 = Eps_1 - (1/3*trEps_1*I);
EpsDev_2 = Eps_2 - (1/3*trEps_2*I);

%produit doublement contracté avec le déviateur
EpsDev1EpsDev1 = doubleProduct(EpsDev_1,EpsDev_1);
EpsDev2EpsDev2 = doubleProduct(EpsDev_2,EpsDev_2);

%% Phase field

d = 0;
g = (1-d)^2;

%% Bourdin

Psi_Bourdin_1 = g*(lambda/2*trEps_1^2 + mu*doubleProduct(Eps_1,Eps_1));
Psi_Bourdin_2 = g*(lambda/2*trEps_2^2 + mu*doubleProduct(Eps_2,Eps_2));

%% Amor

Psi_Amor_1 = g*(K/2*trEpsP_1^2+mu*EpsDev1EpsDev1)+K/2*trEpsM_1^2;
Psi_Amor_2 = g*(K/2*trEpsP_2^2+mu*EpsDev2EpsDev2)+K/2*trEpsM_2^2;

%% Miehe

[Eps_P1, Eps_M1] = SpectralDecomposition(Eps_1);
[Eps_P2, Eps_M2] = SpectralDecomposition(Eps_2);

Psi_Miehe_1 = g*(lambda/2*trEpsP_1^2 + mu*doubleProduct(Eps_P1,Eps_P1)) + ...
                (lambda/2*trEpsM_1^2 + mu*doubleProduct(Eps_M1,Eps_M1));
Psi_Miehe_2 = g*(lambda/2*trEpsP_2^2 + mu*doubleProduct(Eps_P2,Eps_P2)) + ...
                (lambda/2*trEpsM_2^2 + mu*doubleProduct(Eps_M2,Eps_M2));

%% He


%% Resumé

% +---------+--------+--------+
% |         | Zone 1 | Zone 2 |
% +---------+--------+--------+
% | Bourdin |        |        |
% +---------+--------+--------+
% | Amor    |        |        |
% +---------+--------+--------+
% | Miehe   |        |        |
% +---------+--------+--------+


fprintf("\n+---------+------------+------------+");
fprintf("\n| [J/m3]  | Zone 1     | Zone 2     |");
fprintf("\n+---------+------------+------------+");
fprintf("\n| Bourdin |  %.2e  |  %.2e  |",Psi_Bourdin_1, Psi_Bourdin_2);
fprintf("\n+---------+------------+------------+");
fprintf("\n| Amor    |  %.2e  |  %.2e  |",Psi_Amor_1, Psi_Amor_2);
fprintf("\n+---------+------------+------------+");
fprintf("\n| Miehe   |  %.2e  |  %.2e  |",Psi_Miehe_1, Psi_Miehe_2);
fprintf("\n+---------+------------+------------+\n");


%% test

if verif
    % Eps
    test_Eps1 = (2*mu*Eps_1 + lambda*trace(Eps_1)*I)-Sig_1; Test(test_Eps1);
    test_Eps2 = (2*mu*Eps_2 + lambda*trace(Eps_2)*I)-Sig_2; Test(test_Eps2);

    %trace
    test_trEps1 = 3*Sig/E*(2*v-1) - trEps_1; Test(test_trEps1);
    test_trEps2 = Sig/E*(2*v-1) - trEps_2; Test(test_trEps2);
    
    %deviateur
    test_EpsDev1 = [Sig/E*(v+1),0,0;0,-2*Sig/E*(v+1),0;0,0,Sig/E*(v+1)] - EpsDev_1; Test(test_EpsDev1);
    test_EpsDev2 = [Sig/E/3*(v+1),0,0;0,-2*Sig/E/3*(v+1),0;0,0,Sig/E/3*(v+1)] - EpsDev_2; Test(test_EpsDev2);
    
    %produit doublement contracté avec le déviateur
    test_EpsDev1EpsDev1 =  2*(Sig/E*(v+1))^2+(-2*Sig/E*(v+1))^2 - EpsDev1EpsDev1; Test(test_EpsDev1EpsDev1);
    test_EpsDev2EpsDev2 = 2*(Sig/E/3*(v+1))^2+(-2*Sig/E/3*(v+1))^2 - EpsDev2EpsDev2; Test(test_EpsDev2EpsDev2);
    
    % Bourdin
    test_Psi_Bourdin1 = (g*(lambda/2*(3*Sig/E*(2*v-1))^2+mu*9*Sig^2/E^2*(2*v^2+1))) - Psi_Bourdin_1; Test(test_Psi_Bourdin1);
    test_Psi_Bourdin2 = (g*(lambda/2*(Sig/E*(2*v-1))^2+mu*Sig^2/E^2*(2*v^2+1))) - Psi_Bourdin_2; Test(test_Psi_Bourdin2);

    % Amor
    test_Psi_Amor1 = g*mu*(2*(Sig/E*(v+1))^2+(-2*Sig/E*(v+1))^2)+K/2*trEpsM_1^2+ - Psi_Amor_1;
    test_Psi_Amor2 = g*mu*(2*(Sig/E/3*(v+1))^2+(-2*Sig/E/3*(v+1))^2)+K/2*trEpsM_2^2+ - Psi_Amor_2;
    
    % Miehe
    Test(doubleProduct(Eps_P1,Eps_M1));    
    Test(doubleProduct(Eps_P2,Eps_M2))
    Test(Eps_1-(Eps_P1+Eps_M1));
    Test(Eps_2-(Eps_P2+Eps_M2));

    test_Psi_Miehe1 = (g*(mu*2*(9*Sig^2*v^2)/(E^2)) + ((lambda/2*trEps_1^2)+(mu*(9*Sig^2)/(E^2)))) - Psi_Miehe_1;  Test(test_Psi_Miehe1);
    test_Psi_Miehe2 = (g*(mu*2*(Sig^2*v^2)/(E^2)) + ((lambda/2*trEps_2^2)+(mu*(Sig^2)/(E^2)))) - Psi_Miehe_2;  Test(test_Psi_Miehe2);

end

%%

function [Eps_P, Eps_M] = SpectralDecomposition(Eps)
    
    I= eye(3);

    j1 = trace(Eps);
    j2 = (j1^2-trace(Eps^2))/2;
    j3 = det(Eps);

    G = j1^2-3*j2;
    arg = ((2*j1^3) - (9*j1*j2) + (27*j3)) / (2*G^(3/2));
    theta = 1/3*acos(arg);

    e1 = 1/3*j1 + 2/3*G^(1/2)*cos(theta);
    e2 = 1/3*j1 + 2/3*G^(1/2)*cos((2*pi/3)-theta);
    e3 = 1/3*j1 + 2/3*G^(1/2)*cos((2*pi/3)+theta);

    if e1==e2
        e1=1/3*j1 + 1/3*G^(1/2);
        e2=e1;
        e3=1/3*j1 - 2/3*G^(1/2);
        
        M3 = G^(-1/2) * (1/3 * (j1-G^(-1/2))* I-Eps);
        M1M2 = I - M3;

        Eps_P = PositifPart(e1)*M1M2 + PositifPart(e3)*M3;
        Eps_M = NegatifPart(e1)*M1M2 + NegatifPart(e3)*M3;

    elseif e2==e3
        e1=1/3*j1 + 2/3*G^(1/2);
        e2=1/3*j1 - 1/3*G^(1/2);
        e3=e2;

        M1 = G^(-1/2) * (Eps - (1/3*(j1-G^(-1/2))*I));

        M2M3 = I - M1;

        Eps_P = PositifPart(e1)*M1 + PositifPart(e3)*M2M3;
        Eps_M = NegatifPart(e1)*M1 + NegatifPart(e3)*M2M3;
    
    elseif e1==e2 && e2==e3
        e1 = 1/3*j1;

        Eps_P = PositifPart(e1)*I;
        Eps_M = NegatifPart(e1)*I;

    else
        M1 = ((Eps - (e2*I)) / (e1-e2)) * ...
             ((Eps - (e3*I)) / (e1-e3));    
        M2 = ((Eps - (e1*I)) / (e2-e1)) * ...
             ((Eps - (e3*I)) / (e2-e3));
        M3 = ((Eps - (e1*I)) / (e3-e1)) * ...
             ((Eps - (e2*I)) / (e3-e2));
    
        Eps_P = PositifPart(e1)*M1 + PositifPart(e2)*M2 + PositifPart(e3)*M3;
        Eps_M = NegatifPart(e1)*M1 + NegatifPart(e2)*M2 + NegatifPart(e3)*M3;        
    end
       


    

    % [eigenVectors,eigenValues]=eig(Eps_1);

end

function scalar = doubleProduct(A,B)
    scalar = sum(A.*B,"all");
end

function positifPart = PositifPart(a)
    positifPart = (a+abs(a))/2;
end

function negatifPart = NegatifPart(a)
    negatifPart = (a-abs(a))/2;
end

function [] = Test(val)
    tol = 1e-10;
    val=max((abs(val)),[],"all");
    if val>tol
        val
        error('Erreur')
    end
end