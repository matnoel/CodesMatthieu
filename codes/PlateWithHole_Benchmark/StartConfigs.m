% Start all the configs at the same time

clc
clear all
close all

solvers = ["HistoryField"]; %"HistoryField","BoundConstrainedOptim"
splits = ["AnisotropicMiehe", "AnisotropicHe"]; % "Isotropic", "AnisotropicAmor", "AnisotropicMiehe", "AnisotropicHe"
regularizations = ["AT1", "AT2"]; % "AT1", "AT2"

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

