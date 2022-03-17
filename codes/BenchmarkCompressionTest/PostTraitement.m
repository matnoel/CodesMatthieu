
clc
clearvars
close all

% Pour tout les modèles il faut :
% - Afficher l'endommagement pour : ud = 18.5 µm et ud = 24.6µm
% - Afficher courbe déplacement
% - Afficher les temps de résolution

test = true;
% test = false;

solvers = ["HistoryField","BoundConstrainedOptim"]; %"HistoryField","BoundConstrainedOptim"
splits = ["AnisotropicMiehe", "AnisotropicHe"]; % "Isotropic", "AnisotropicAmor", "AnisotropicMiehe", "AnisotropicHe"
regularizations = ["AT1","AT2"]; % "AT1", "AT2"

% Fichier dans lequel on va récupérer les données
foldername = 'CompressionTest_PlateWithHole';

% Dossier ou on sauvegarde le post traitement
filename = 'PostTraitement';
pathnamePostTraitement = BuiltPathnameResult(test,foldername,filename);

if ~exist(pathnamePostTraitement,'dir')
    mkdir(pathnamePostTraitement);
end

%%

legends = {};
i=0;

% Pour chaque solveur
for s=1:length(solvers)    
    solver = solvers(s);
    
    % Pour chaque modèle Phase field
    for sp=1:length(splits)
        split = splits(sp);
        
        % Pour chaque modèle de régularisation
        for r=1:length(regularizations)
            regularization = regularizations(r);

            filename = append(solver,'_',split,'_', regularization);
            pathname = BuiltPathnameResult(test,foldername,filename);

            if ~exist(pathname,'dir')
                break
            end

            i=i+1;

            % Charge le problème
            load(fullfile(pathname,'solution.mat'),'phaseFieldSolution');
            
            if iscell(phaseFieldSolution)
                phaseFieldSolution = phaseFieldSolution{1};
            end

            PFM = phaseFieldSolution.PFM;
            S = PFM.S;
            S_phase = PFM.S_phase;

            udt = phaseFieldSolution.udt;
            ft = phaseFieldSolution.ft;
            dt = phaseFieldSolution.dt;
            resolutionTime = phaseFieldSolution.resolutionTime;

            temps = GetTime(sum(resolutionTime));

            fprintf("\n"+filename+" - "+temps+"\n");

            %% Affichage 
            
            % Force déplacement
            figure(1)
            plot(udt*1e6,-ft*1e-6,'LineWidth',1)
            hold on
            legends{i} = append(split,' ', regularization);
            
            % Endommagement
            plotSolution(S_phase,dt{length(dt)});
            filename = append('damage ',filename,'.png');
            saveas(gcf, fullfile(pathnamePostTraitement, filename))
        end
        
    end
    
    figure(1)
    grid on
    title(solver,'interpreter','Latex','fontsize',15)
    xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
    ylabel("Load in kN/mm",'interpreter','Latex','fontsize',15)  
    legend(legends,'Location','southeast','interpreter','Latex','fontsize',13)
    file = append('displacement ',solver,'.png');
    saveas(gcf, fullfile(pathnamePostTraitement, file))
    hold off

    i=0;

    close all
end



fprintf('\n plot saved \n')