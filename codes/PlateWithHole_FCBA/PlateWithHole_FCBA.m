%% Benchmark Compression test - Plate with hole %%
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]

% clc
clearvars
close all

%% Options

postTraitement = false;

% test = true;
test = false;

if postTraitement
    setPFM = false;
    solve = false;
    plotResults = true;
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

foldername = 'PlateWithHole_FCBA';

solver = 'HistoryField'; %'HistoryField','BoundConstrainedOptim'
split = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
regularization = 'AT1'; % 'AT1', 'AT2'

depSurCercle = true;

filename = append(solver,'_', split,'_',regularization);

if depSurCercle
    filename = append(filename,'_cercle');
end

pathname = BuiltPathnameResult(test,foldername,filename);
    
if ~exist(pathname,'dir')
    pathname
    mkdir(pathname);
end

%% Set Phase Field Model

if setPFM

    %% Datas Geo
    
    H = 120e-3; %[m]
    L = 90e-3;
    ep = 20e-3;   
    
    holeDiameter = 10e-3;

    h=35e-3;
    
    %% Datas Phase Field
    
    gc = 2.7/3/2.68/1.85;   % Taux de libération d'énergie critique (ou ténacité à la rupture) [N/m]
%     l_0 = 1.2e-3; % Paramètre de régularisation (largeur de la fissure)
    l_0 = L/100; % Paramètre de régularisation (largeur de la fissure)
    
    %% Création du modèle et maillage  

    clC = l_0/2;
    clD = clC*10;    

    if test
        clD = 2*l_0;
        clC = l_0;
    end

    domain = DOMAIN(2,[-L/2,0],[L/2,H]);
    circle = CIRCLE(0,H-h,holeDiameter/2);

    optimesh = true;
    
    model = gmshplateWithHole(domain,circle,clD,clC,split,optimesh,pathname);
        
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
    E = 10e9;
    E=E/160*1.8;
%     E= 11.580*1e9;
    NU = 0.3;
            
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',split);
    
    mat = setnumber(mat,1);
    
    S = setoption(model,'CONT');
    S = setmaterial(S,mat);
    
    %% Conditions en déplacement
    
    BLower = LIGNE([-L/2, 0],[L/2, 0]);
    BUpper = LIGNE([-L/2, H],[L/2, H]);
    BLeft = LIGNE([-L/2, 0], [-L/2, H]); 
    BRight = LIGNE([L/2, 0], [L/2, H]);   
    
    P0 = getvertices(domain);
    P0 = POINT(P0{1});

    S = final(S);

    % Get the nodes on the egdes
    [~,noeudsBLower] = intersect(S, BLower);
    [~,noeudsBUpper] = intersect(S, BUpper);
    [~,noeudsBLeft] = intersect(S, BLeft);
    [~,noeudsBRight] = intersect(S, BRight);
    noeudsDuBord = [noeudsBLeft; noeudsBRight; noeudsBLower; noeudsBUpper];

    % Nodes of the circle
    [~,nodesBCircle] = intersect(S,circle);  
    coord = getcoord(S.node(nodesBCircle));
    
    % nodes of the lower circle
    if depSurCercle
        loadNodes = nodesBCircle(find(coord(:,2)-(H-h)<1e-6));
    else
        [~, loadNodes] = intersect(S,BUpper);
    end

    phaseFieldModel.DirichletBoundaryConditions{1} = {BLower, 'UY', 0};
    phaseFieldModel.DirichletBoundaryConditions{2} = {P0, 'UX', 0};
    % apply the boundary conditions
    S = ApplyDirichletBoundaryConditions(S, phaseFieldModel.DirichletBoundaryConditions);
    
    phaseFieldModel.S = S;

    % the last conditions must contain loadnodes and directions for phase
    % field loading
    phaseFieldModel.DirichletBoundaryConditions{3} = {loadNodes, 'UY', 0};

%     %plot mesh
%     plotModel(model,'Color','k','FaceColor','c','FaceAlpha',1,'legend',false);
%     hold on3
%     coordMesh = getcoord(S.node);
%     plot(coordMesh(loadNodes,1), coordMesh(loadNodes,2), 'rx', LineWidth=1.5)

    %% Load
            
%     inc0 = 8e-8;   %[m] step d<=0.6
%     inc1 = 2e-8;   % step d>0.6
%     umax = 25e-6; % 25 µm
%     
%     if test     
%         inc0 = 16e-8;
%         inc1 = 4e-8;
%     end    

    umax = 3e-3; % 3mm
    Ninc = 400;
    inc0 = umax/Ninc; %[m] step d<=0.6
    inc1 = inc0/3; % step d>0.6

    if test
        inc0=inc0*2;
        inc1=inc1*2;                
    end
    
    displacement = struct('inc0',inc0,'inc1',inc1,'umax',umax,'dthreshold',0.6);

    save(fullfile(pathname,'problem.mat'), 'phaseFieldModel','displacement');
    
else
    load(fullfile(pathname,'problem.mat'),'phaseFieldModel','displacement');
%     fprintf(phaseFieldModel.resume)
end

%% Resolution

if solve
%     phaseFieldSolution = PhaseFieldTresholdSimulation(PFM, loadNodes,inc1, inc0, umax, 0.6, display);
    signe = -1;
    phaseFieldSolution = PhaseFieldSimulation(phaseFieldModel, displacement,"mm",signe, display,noeudsDuBord);
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
