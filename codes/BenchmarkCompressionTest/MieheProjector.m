% ici j'essaye de construire le projecteur

clc, close all; clear all

Eps=[1,-3;-3,2]

[eigenVectors,eigenValues]=eig(Eps);

% invariants 1 et 2
I1=trace(Eps);  
I2=det(Eps);

%valeurs propre 1 et 2
l1 = (I1-sqrt(I1^2-4*I2))/2
l2 = (I1+sqrt(I1^2-4*I2))/2

norm1=sqrt((-Eps(1,2)/(Eps(1,1)-l1))^2+1^2);
norm2=sqrt((-Eps(1,2)/(Eps(1,1)-l2))^2+1^2);

%vecteur prorres 1 et 2
n1=1/norm1*[-Eps(1,2)/(Eps(1,1)-l1);1] %[x1,y1]
x1=n1(1,1);
y1=n1(2,1);

n2=1/norm2*[-Eps(1,2)/(Eps(1,1)-l2);1] %[x2,y2]
x2=n2(1,1);
y2=n2(2,1);

%test si la décomposition marche
test1=Eps*n1-n1*l1;
test2=Eps*n2-n2*l2;

% test orthognolnalité
test3=sum(n1.*n2);

%base des vecteurs propres
M1=[x1*x1,x1*y1;y1*x1,y1*y1]
M2=[x2*x2,x2*y2;y2*x2,y2*y2]

% test orthognolnalité
test4=sum(M1.*M2,"all");

%test si la base est identique
M1_verif=(Eps-l2*eye(2))/(l1-l2);
M2_verif=eye(2)-M1_verif;

test_M1=M1-M1_verif;
test_M2=M2-M2_verif;


