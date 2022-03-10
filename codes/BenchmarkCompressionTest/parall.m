
clc
clearvars
close all

solvers = ["BoundConstrainedOptim"]; %"HistoryField","BoundConstrainedOptim"
PFmodels = ["AnisotropicMiehe"]; % "Isotropic", "AnisotropicAmor", "AnisotropicMiehe", "AnisotropicHe"
alphas = ["AT1","AT2"]; % "AT1", "AT2"

listPar=[];

Np=0;
for s=1:length(solvers)
    for pf=1:length(PFmodels)
        for a=1:length(alphas)
            Np = Np+1;
            listPar{Np} = [solvers(s), PFmodels(pf), alphas(a)]            
        end
    end    
end

Np

delete(gcp('nocreate'))
parpool(Np)

parfor p=1:Np
    
    solver = listPar{p}{1};
    PFmodel = listPar{p}{2};
    alpha = listPar{p}{3};

    BenchmarkCompressionTest(solver, PFmodel, alpha);
end

delete(gcp('nocreate'))