clc
clearvars
close all

display = true;
makeMovie = false;
plot = true;
saveParaview = true;

symmetry = 'Isotropic'; % 'Isotropic' or 'Anisotropic'. Material symmetry
% option = 'DEFO'; % plane strain 
option = 'CONT'; % plane stress

filename = append('PoutreFlexion');

pathname = fullfile(getfemobjectoptions('path'),'MYCODE','results','codeMatthieu','PoutreFlexion',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

%% Datas

timTot = tic;

Dim = 2;

P = 800;

Length = 120; %[mm]
height = 13;
ep = height;

%% Création du modèle et maillage

nby = 10;
cl = height/nby; %longeur de maille
nby = ceil(Length/cl);

tMesh = tic;

domain = DOMAIN(Dim,[0,0],[Length,height]);

elemtype = 'TRI3'; %'QUA4', 'TRI6', 'QUA4', 'QUA8'
model = build_model(domain,'cl',cl,'elemtype',elemtype,'option',option,'recombine','filename',fullfile(pathname,'gmsh_domain'));

% % nbelem = repmat(10,1,Dim);
% nbelem = [nbx; nby];
% model = build_model(domain,'nbelem',nbelem,'elemtype',elemtype,'option',option);

fprintf('Ne = %d, Nn = %d \n', model.nbelem, model.nbnode)

timeMesh = toc(tMesh);
fprintf('Meshing: elapsed time = %f s\n',timeMesh);

% plotModel(model,'Color','k','FaceColor','c','FaceAlpha',1,'legend',false);

%% Matériau

E = 210000.0;
NU = 0.3;
RHO = 1;

mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep);
mat = setnumber(mat,1);
model = setmaterial(model, mat);


%% Boundary conditions

tAssemble = tic;

B0 = LIGNE([0.0,0],[0,height]);
BL = LIGNE([Length,0.0],[Length,height]);
    
model = final(model);
model = addcl(model,B0);

[~, loadNodes] = intersect(model,BL);

charge = P;
chargeSurfacique = charge/(height);
% charge = charge/length(loadNodes);

b = surfload(model ,BL,'FY',-chargeSurfacique);

Kglob = calc_rigi(model);

timeAssemble = toc(tAssemble);
fprintf('Assembling: elapsed time = %f s\n',timeAssemble);

%% Solve

tSolve = tic;

u = Kglob\b;
min_uy = min(u)

% R = chol(Kglob);
% u = R\(R'\b);

timeSolve = toc(tSolve);
fprintf('Solving: elapsed time = %f s\n',timeSolve);

tPostProcess = tic;
Wdef = 1/2*u'*Kglob*u
% Wdef = Wdef(1)
e = calc_epsilon(model,u);

% Contraintes
stress = double(mean(calc_sigma(model,u),4));

Sxx = reshape(stress(1,1,:),[getnbelem(model),1]);
Syy = reshape(stress(2,1,:),[getnbelem(model),1]);
Sxy = reshape(stress(3,1,:),[getnbelem(model),1]);

Svm = sqrt(Sxx.^2+Syy.^2-(Sxx.*Syy)+3*Sxy.^2);

Svm_max = max(Svm)/13

timePostProcess = toc(tPostProcess);
fprintf('Post-processing: elapsed time = %f s\n',timePostProcess);

timTot = toc(timTot);
fprintf('Temps total = %f s\n',timTot);

%% Post Process

plotDomain(domain,'legend',false);

plotModel(model,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);

[hD,legD] = plotBoundaryConditions(model,'legend',false);
[hN,legN] = vectorplot(model,'F',b,'r','LineWidth',1);
            

ampl = getsize(model)/max(abs(u))/5;

for i=1:Dim
    plotSolution(model,u,'displ',i,'ampl',ampl);
end

% for i=1:(Dim*(Dim+1)/2)
%     plotSolution(model,u,'epsilon',i,'ampl',ampl);
%     plotSolution(model,u,'sigma',i,'ampl',ampl);    
% end

% plotSolution(model,u,'sigma','mises','ampl',ampl);

