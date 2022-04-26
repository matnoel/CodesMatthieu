
% clc
clearvars
close all

%% Options

postTraitement = false;

test = true;
% test = false;

if postTraitement
    setPFM = false;
    solve = false;
    plotResults = false;
    saveParaview = true;
    makeMovie = false;
else
    setPFM = true;
    solve = true;
    plotResults = true;
    saveParaview = true;
    makeMovie = false;
end

display = true;

%% Configs

foldername = 'ShearTest';

solver = 'HistoryField'; %'HistoryField','BoundConstrainedOptim'
split = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
regularization = 'AT2'; % 'AT1', 'AT2'

filename = append(solver,'_', split,'_',regularization);

pathname = BuiltPathnameResult(test,foldername,filename);
    
if ~exist(pathname,'dir')
    pathname
    mkdir(pathname);
end

%% Set Phase Field Model

if setPFM

    %% Datas Geo
    
    DIM = 2;
    
    L = 1e-3;
    a = L/2;
    ep = 1;

    D = DOMAIN(2,[0.0,0.0],[L,L]);
    C = LIGNE([0.0,L/2],[a,L/2]);
    
    %% Datas Phase Field
    
    gc = 2.7e3;   % Taux de libération d'énergie critique (ou ténacité à la rupture) [N/m]
    l_0 = 1e-5; % Paramètre de régularisation (largeur de la fissure)
    
    %% Création du modèle et maillage

    clD = l_0/2; %clD = 2.5e-5;
    clC = clD; % clC = 2.5e-6;
    
    if test
        clD = 1e-5;
        clC = 1e-5;
    end

    gmshfile = fullfile(pathname,'gmsh_domain_single_edge_crack');
    model = gmshdomainwithedgecrack(D,C,clD,clC,gmshfile);

%     S_phase = gmshdomainwithedgesmearedcrack(D,C,c,clD,clC,fullfile(pathname,'gmsh_domain_single_edge_crack'));
        
    %% Modele a gradient d'endommagement
    
    g = @(d) (1-d).^2;  % Fonction degradation energetique    
    
    phaseFieldModel = PhaseFieldModel(solver, split, regularization, gc, l_0);
    
    % Materiau
    mat_phase = FOUR_ISOT('k',phaseFieldModel.k,'r',phaseFieldModel.r(0));
    mat_phase = setnumber(mat_phase,1);

    % Attribu le matériau
    S_phase = setmaterial(model,mat_phase);
    
    % Finalise et applique Cl
%     S_phase = final(S_phase);
    S_phase = final(S_phase);
    S_phase = addcl(S_phase,C,'T',1);

    
    
    phaseFieldModel.S_phase = S_phase;
    
    d = calc_init_dirichlet(S_phase);    
    
    %% Modele deformation elastique
    
    RHO = 1;    % Densité kg/m3
    k = 1e-10;  % Petite rigidité résiduelle artificielle
    E = 210e9;
    NU = 0.3;
            
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',split);
    
    mat = setnumber(mat,1);
    
    S = setoption(model,'DEFO'); % 'CONT'
    S = setmaterial(S,mat);
    
    %% Conditions en déplacement
    
    BLower = LIGNE([0, 0],[L, 0]);
    BUpper = LIGNE([0, L],[L, L]);
    BLeft = LIGNE([0, 0], [0, L]); 
    BRight = LIGNE([L, 0], [L, L]);

    S = final(S);

    [~, loadNodes] = intersect(S,BUpper);

    phaseFieldModel.DirichletBoundaryConditions{1} = {BLower, {'UX','UY'}, [0;0]};
    phaseFieldModel.DirichletBoundaryConditions{2} = {BLeft, 'UY', 0};
    phaseFieldModel.DirichletBoundaryConditions{3} = {BRight, 'UY', 0};
%     phaseFieldModel.DirichletBoundaryConditions{4} = {BUpper, 'UY', 0};
    
    % apply the boundary conditions
    S = ApplyDirichletBoundaryConditions(S, phaseFieldModel.DirichletBoundaryConditions);
    
    phaseFieldModel.S = S;

    % the last conditions must contain loadnodes and directions for phase
    % field loading
    phaseFieldModel.DirichletBoundaryConditions{4} = {BUpper, 'UX', 0};

    %% Load
            
    dt = 1e-8;
    % nt = 1500;
    nt = 2000;
    if test
        dt = 5e-8;
        % nt = 300;
        nt = 400;
    end
    displacement = linspace(dt,nt*dt,nt);

    displacement = TIMEMODEL(displacement);
    
    save(fullfile(pathname,'problem.mat'), 'phaseFieldModel','displacement');
    
else
    load(fullfile(pathname,'problem.mat'),'phaseFieldModel','displacement');
%     fprintf(phaseFieldModel.resume)
end

%% Resolution

if solve
    phaseFieldSolution = PhaseFieldSimulation(phaseFieldModel, displacement,"mm",1, display,[]);
    save(fullfile(pathname,'solution.mat'),'phaseFieldSolution');    
else
    load(fullfile(pathname,'solution.mat'),'phaseFieldSolution');    
end

%% Save solutions

if plotResults
    phaseFieldSolution.PlotResults(pathname,"mm");
end

if saveParaview   
    phaseFieldSolution.SaveParaview(pathname);
end

if makeMovie
    phaseFieldSolution.MakeMovie(pathname)
end
