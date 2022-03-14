%% Benchmark Compression test - Plate with hole %%
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]

clc
clearvars
close all

%% Options

parallel = false;

postTraitement = false;

% test = true;
test = false;

if postTraitement
    setPFM = false;
    solve = false;
    plotResults = true;
    saveParaview = false;
else
    setPFM = true;
    solve = true;
    plotResults = true;
    saveParaview = false;
end


display = true;

%% Configs

foldername = 'CompressionTest_PlateWithHole';

solvers = ["HistoryField"]; %"HistoryField","BoundConstrainedOptim"
splits = ["AnisotropicMiehe"]; % "Isotropic", "AnisotropicAmor", "AnisotropicMiehe", "AnisotropicHe"
regularizations = ["AT1"]; % "AT1", "AT2"

configs={};
nbConfig=0;

for s=1:length(solvers)
    for sp=1:length(splits)
        for r=1:length(regularizations)
            nbConfig = nbConfig+1;
            configs{nbConfig} = [solvers(s), splits(sp), regularizations(r)];

            % Construct folder
            filenames{nbConfig} = append(solvers(s),'_', splits(sp),'_',regularizations(r));

            pathnames{nbConfig} = BuiltPathnameResult(test,foldername,filenames{nbConfig});
    
            if ~exist(pathnames{nbConfig},'dir')
                mkdir(pathnames{nbConfig});
            end

        end
    end    
end
   

%% Set Phase Field Model

% for each configuration we construct the geometry, the phase field model
% and we save in the good folder
% or we load all the models

if setPFM
    for config=1:nbConfig
        
        pathname = pathnames{config};
        solver = configs{config}{1};
        split = configs{config}{2};
        regularization = configs{config}{3};
    
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
        
        PFM = PhaseFieldModel(solver, split, regularization, gc, l_0);
        
        % Materiau
        mat_phase = FOUR_ISOT('k',PFM.k,'r',PFM.r(0));
        mat_phase = setnumber(mat_phase,1);
    
        % Attribu le matériau
        S_phase = setmaterial(model,mat_phase);
        
        % ici pas besoin de de renseigner des conditions limite car pas de fissure preexitante
        S_phase = final(S_phase);
        
        PFM.S_phase = S_phase;
        
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
        
        [~, loadNodes] = intersect(S,BUpper);
        
        % Get the nodes on the egdes
        [~,noeudsBLower] = intersect(S, BLower);
        [~,noeudsBUpper] = intersect(S, BUpper);
        [~,noeudsBLeft] = intersect(S, BLeft);
        [~,noeudsBRight] = intersect(S, BRight);
        noeudsDuBord = [noeudsBLeft; noeudsBRight; noeudsBLower; noeudsBUpper];
        
        S = final(S);
        PFM.DirichletBoundaryConditions{1} = {BLower, 'UY', 0};
        PFM.DirichletBoundaryConditions{2} = {[0,0], 'UX', 0};
        
        S = ApplyDirichletBoundaryConditions(S, PFM.DirichletBoundaryConditions);
    
        PFM.S = S;
        PhaseFieldModels{config} = PFM;    
        save(fullfile(pathname,'problem.mat'), 'PFM', 'loadNodes');
    end
else
    for config=1:nbConfig
        load(fullfile(pathnames{config},'problem.mat'));
        PhaseFieldModels{config} = PFM;        
    end    
end

%% Load
        
bigInc = 8e-8;   %[m] step d<=0.6
smallInc = 2e-8;   % step d>0.6
maxDep = 25e-6; % 25 µm

if test     
    bigInc = 16e-8;
    smallInc = 4e-8;
end

%% Resolution

if solve   
    
    phaseFieldSolutions = {};

    if parallel
        delete(gcp('nocreate'))
        parpool(nbConfig)
        
        parfor config=1:nbConfig 
            PFM=PhaseFieldModels{config};
            phaseFieldSolutions{config} = PhaseFieldTresholdSimulation(PFM, loadNodes,smallInc, bigInc, maxDep, 0.6, display);
        end

        delete(gcp('nocreate'))
    else        
        for config=1:nbConfig
            PFM=PhaseFieldModels{config};
            phaseFieldSolutions{config} = PhaseFieldTresholdSimulation(PFM, loadNodes,smallInc, bigInc, maxDep, 0.6, display);            
        end
    end

    % Save
    for config=1:nbConfig
        phaseFieldSolution = phaseFieldSolutions{config};
        save(fullfile(pathnames{config},'solution.mat'),'phaseFieldSolution');
    end

    
else
    for config=1:nbConfig
        load(fullfile(pathnames{config},'solution.mat'),'phaseFieldSolution');
        phaseFieldSolutions{config} = phaseFieldSolution;
    end
end

%% Save solutions

if plotResults || saveParaview
    for config=1:nbConfig       
        
        pathname = pathnames{config};
        phaseFieldSolution = phaseFieldSolutions{config};
        PFM = PhaseFieldModels{config};

        if plotResults
            phaseFieldSolution.PlotResults(pathname)
        end
    
%         if saveParaview   
%             SaveParaview(ut, dt, PFM);
%         end
    end
end
