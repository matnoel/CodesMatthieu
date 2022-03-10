
clc
clearvars
close all

solvers = ["BoundConstrainedOptim"]; %"HistoryField","BoundConstrainedOptim"
splits = ["AnisotropicMiehe"]; % "Isotropic", "AnisotropicAmor", "AnisotropicMiehe", "AnisotropicHe"
regularizations = ["AT1","AT2"]; % "AT1", "AT2"

listPar={};

Np=0;
for s=1:length(solvers)
    for sp=1:length(splits)
        for r=1:length(regularizations)
            Np = Np+1;
            listPar{Np} = [solvers(s), splits(sp), regularizations(r)];
        end
    end    
end

Np

delete(gcp('nocreate'))
parpool(Np)

parfor p=1:Np
    
    solver = listPar{p}{1};
    split = listPar{p}{2};
    regularization = listPar{p}{3};

    BenchmarkCompressionTest(solver, split, regularization);
end

delete(gcp('nocreate'))