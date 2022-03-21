%% Benchmark Compression test - Plate with hole %%
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]

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
else
    setPFM = true;
    solve = true;
    plotResults = true;
    saveParaview = true;
end

display = true;

%% Configs

foldername = 'CompressionTest_PlateWithHole';

solver = 'HistoryField'; %'HistoryField','BoundConstrainedOptim'
split = 'AnisotropicHe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
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
    % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
    
    DIM = 2;
    
    h = 30e-3; %[m]
    L = 15e-3;
    ep = 1;   
    
    holeDiameter = 12e-3/2;
    
    %% Datas Phase Field
    
    gc = 1.4;   % Taux de libération d'énergie critique (ou ténacité à la rupture) [N/m]
    l_0 = 0.12e-3; % Paramètre de régularisation (largeur de la fissure)
    
    %% Création du modèle et maillage  

    clD = l_0/2;
    clC = clD;

    if test
        clD = 0.25e-3;
        clC = 0.12e-3;
    end

    domain = DOMAIN(2,[0,0],[L,h]);
    circle = CIRCLE(L/2,h/2,holeDiameter/2);
    gmshfile = fullfile(pathname,'plateWithHole');
    model = gmshdomainwithhole(domain,circle,clD,clC,gmshfile);
        
    %% Modele a gradient d'endommagement
    
    g = @(d) (1-d).^2;  % Fonction degradation energetique    
    
    phaseFieldModel = PhaseFieldModel(solver, split, regularization, gc, l_0);
    
    % Materiau
    mat_phase = FOUR_ISOT('k',phaseFieldModel.k,'r',phaseFieldModel.r(0));
    mat_phase = setnumber(mat_phase,1);

    % Attribu le matériau
    S_phase = setmaterial(model,mat_phase);
    
    % ici pas besoin de de renseigner des conditions limite car pas de fissure preexitante
    S_phase = final(S_phase);
    
    phaseFieldModel.S_phase = S_phase;
    
    d = calc_init_dirichlet(S_phase);    
    
    %% Modele deformation elastique
    
    RHO = 1;    % Densité kg/m3
    k = 1e-10;  % Petite rigidité résiduelle artificielle
    E = 12e9;
    NU = 0.3;
            
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',split);
    
    mat = setnumber(mat,1);
    
    S = setoption(model,'DEFO');
    S = setmaterial(S,mat);
    
    %% Conditions en déplacement
    
    BLower = LIGNE([0, 0],[L, 0]);
    BUpper = LIGNE([0, h],[L, h]);
    BLeft = LIGNE([0, 0], [0, h]); 
    BRight = LIGNE([L, 0], [L, h]);   
    
    P0 = getvertices(domain);
    P0 = POINT(P0{1});

    S = final(S);

    [~, loadNodes] = intersect(S,BUpper);
    
    % Get the nodes on the egdes
    [~,noeudsBLower] = intersect(S, BLower);
    [~,noeudsBUpper] = intersect(S, BUpper);
    [~,noeudsBLeft] = intersect(S, BLeft);
    [~,noeudsBRight] = intersect(S, BRight);
    noeudsDuBord = [noeudsBLeft; noeudsBRight; noeudsBLower; noeudsBUpper];   

    phaseFieldModel.DirichletBoundaryConditions{1} = {BLower, 'UY', 0};
%     phaseFieldModel.DirichletBoundaryConditions{2} = {[0,0], 'UX', 0};
    phaseFieldModel.DirichletBoundaryConditions{2} = {P0, 'UX', 0};
    % apply the boundary conditions
    S = ApplyDirichletBoundaryConditions(S, phaseFieldModel.DirichletBoundaryConditions);
    
    phaseFieldModel.S = S;

    % the last conditions must contain loadnodes and directions for phase
    % field loading
    phaseFieldModel.DirichletBoundaryConditions{3} = {loadNodes, 'UY', 0};

    %% Load
            
    inc0 = 8e-8;   %[m] step d<=0.6
    inc1 = 2e-8;   % step d>0.6
    umax = 25e-6; % 25 µm
    
    if test     
        inc0 = 16e-8;
        inc1 = 4e-8;
    end
    
    displacement = struct('inc0',inc0,'inc1',inc1,'umax',umax,'dthreshold',0.6);

    save(fullfile(pathname,'problem.mat'), 'phaseFieldModel','displacement');
    
else
    load(fullfile(pathname,'problem.mat'),'phaseFieldModel','displacement');
end

%% Resolution

if solve
%     phaseFieldSolution = PhaseFieldTresholdSimulation(PFM, loadNodes,inc1, inc0, umax, 0.6, display);
    phaseFieldSolution = PhaseFieldSimulation(phaseFieldModel, displacement, display);
    save(fullfile(pathname,'solution.mat'),'phaseFieldSolution');
else
    load(fullfile(pathname,'solution.mat'),'phaseFieldSolution');
end

%% Save solutions

if plotResults || saveParaview
    
    if plotResults
        phaseFieldSolution.PlotResults(pathname);
    end

    if saveParaview   
        phaseFieldSolution.SaveParaview(pathname);
    end

end

