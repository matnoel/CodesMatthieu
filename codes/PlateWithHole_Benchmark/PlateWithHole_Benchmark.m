%% Benchmark Compression test - Plate with hole %%
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]

% clc
clearvars
close all

%% Options

postTraitement = true;

% test = true;
test = false;

if postTraitement
    setPFM = false;
    solve = false;
    plotResults = true;
    saveParaview = false;
    makeMovie = false;
else
    setPFM = true;
    solve = true;
    plotResults = true;
    saveParaview = false;
    makeMovie = false;
end

display = true;

%% Configs

foldername = 'PlateWithHole_Benchmark';

solver = 'HistoryField'; %'HistoryField','BoundConstrainedOptim'
split = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
regularization = 'AT1'; % 'AT1', 'AT2'

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
    
    coef = 1e-3;
%     coef=1;


    h = 30*coef; %[m]
    L = 15*coef;
    ep = 1;   
    
    holeDiameter = 6 * coef;
    
    %% Datas Phase Field
    
    gc = 1.4;   % Taux de libération d'énergie critique (ou ténacité à la rupture) [N/m]
%     gc = 1.4/1000;   % [N/mm]
    l_0 = 0.12*coef; % Paramètre de régularisation (largeur de la fissure)
    
    %% Création du modèle et maillage  

    clD = l_0/2;
    clC = clD;

    if test
        clD = 0.25*coef;
        clC = 0.12*coef;
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
    E = 12e9;   %Pa
%     E = 12e9*1e-6;   %MPa
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
%     fprintf(phaseFieldModel.resume)
end

%% Resolution

if solve
%     phaseFieldSolution = PhaseFieldTresholdSimulation(PFM, loadNodes,inc1, inc0, umax, 0.6, display);
    phaseFieldSolution = PhaseFieldSimulation(phaseFieldModel, displacement,"µm",-1, display);
    save(fullfile(pathname,'solution.mat'),'phaseFieldSolution');    
else
    load(fullfile(pathname,'solution.mat'),'phaseFieldSolution');    
end

%% Save solutions

if plotResults
    unite = "µm";
    phaseFieldSolution.PlotResults(pathname, unite);
end

if saveParaview   
    phaseFieldSolution.SaveParaview(pathname);
end

if makeMovie
    phaseFieldSolution.MakeMovie(pathname)
end

%% Look at the stresses and strains and energy near the hole
% 
% clc
% 
% [~, nodesOnCircle] = intersect(S,circle);
% coord = getcoord(S.node(nodesOnCircle));
% 
% dec = 1e-5;
% 
% nodes1 = find(coord(:,1) >= L/2 + holeDiameter/2 - dec);
% nodes2 = find(coord(:,2) >= L + holeDiameter/2 - dec);
% 
% plotModel(model,'Color','k','FaceColor','c','FaceAlpha',1,'legend',false);
% hold on 
% scatter(coord(nodes1,1),coord(nodes1,2),'bx','LineWidth',1.5)
% scatter(coord(nodes2,1),coord(nodes2,2),'rx','LineWidth',1.5)
% 
% elems = getgroupelem(S);
% connect = getconnec(elems{1});
% 
% elems1 = GetElements(connect, nodes1);
% elems2 = GetElements(connect, nodes2);
% 
% % for i=1:length(phaseFieldSolution.dt)
% % 
% %     ui = phaseFieldSolution.ut{i};
% %     di = phaseFieldSolution.dt{i};
% %     
% %     % update mat
% %     mats = MATERIALS(S);
% %     for m=1:length(mats)
% %         mats{m} = setparam(mats{m},'d',di);
% %         mats{m} = setparam(mats{m},'u',ui);
% %     end
% %     S = actualisematerials(S,mats);
% %     
% %     % Stress
% %     stress = double(mean(calc_sigma(S,ui),4));
% %     
% %     Sxx = reshape(stress(1,1,:),[getnbelem(S),1]);
% %     Syy = reshape(stress(2,1,:),[getnbelem(S),1]);
% %     Sxy = reshape(stress(3,1,:),[getnbelem(S),1]);
% % 
% %     Svm = sqrt(Sxx.^2+Syy.^2-(Sxx.*Syy)+3*Sxy.^2);
% % 
% %     Sxx_t_e(i,1) = mean(Sxx(elems1)); Sxx_t_e(i,2) = mean(Sxx(elems2));
% %     Syy_t_e(i,1) = mean(Syy(elems1)); Syy_t_e(i,2) = mean(Syy(elems2));
% %     Sxy_t_e(i,1) = mean(Sxy(elems1)); Sxy_t_e(i,2) = mean(Sxy(elems2));
% % 
% % end
% % 
% % fprintf("")
% 
% %% 
% 
% function elems = GetElements(connect, nodes)
%     % Get the elements who are using the nodes
%     
%     [Lia, ~] = ismember(connect, nodes);
%     [row,~] = find(Lia>0);
%     elems = unique(row);
% end