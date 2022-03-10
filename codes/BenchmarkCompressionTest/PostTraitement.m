
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

% Dossier ou on sauvegarde le post traitement
PostTraitement_filename = append('BenchmarkCompressionTest_PostTraitement');

if test
    PostTraitement_filename = append(PostTraitement_filename,'_Test');
    PostTraitement_pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest','Test',PostTraitement_filename);
else
    PostTraitement_pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest',PostTraitement_filename);
end

if ~exist(PostTraitement_pathname,'dir')
    mkdir(PostTraitement_pathname);
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

            name = append(solver,' ',split,' ', regularization);

            filename = append('BenchmarkCompressionTest_',solver,'_', split,'_',regularization);

            if test
                filename = append(filename,'_Test');
                pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest','Test',filename);
            else
                pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest',filename);        
            end

            if ~exist(pathname,'dir')
                break
            end

            i=i+1;

            % Charge le problème
            load(fullfile(pathname,'problem.mat'),'PFM');    
            S = PFM.S;
            S_phase =PFM.S_phase;

            % Charge la solution (que ce qui nous intéresse)
            load(fullfile(pathname,'solution.mat'),'ud_t','dt','ut','ft','resolutionTime'); 
            
            temps = GetTime(sum(resolutionTime));

            fprintf("\n"+name+" - "+temps+"\n");

            %% Affichage 
            
            % Force déplacement
            figure(1)
            plot(ud_t*1e6,-ft*1e-6,'LineWidth',1)
            hold on
            legends{i} = append(split,' ', regularization);
            
            % Endommagement
            plotSolution(S_phase,dt{length(dt)});
            name = append('damage ',name,'.png');
            saveas(gcf, fullfile(PostTraitement_pathname, name))

        end
    end
    
    figure(1)
    grid on
    title(solver,'interpreter','Latex','fontsize',15)
    xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
    ylabel("Load in kN/mm",'interpreter','Latex','fontsize',15)  
    legend(legends,'Location','southeast','interpreter','Latex','fontsize',13)
    file = append('displacement ',solver,'.png');
    saveas(gcf, fullfile(PostTraitement_pathname, file))
    hold off
end



fprintf('\n plot saved \n')