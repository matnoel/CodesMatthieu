function [] = BenchmarkCompressionTest(solver, PFmodel, alpha)

%% Benchmark Compression test %%
% [Nguyen, Yvonnet, Waldmann, He, 2020, IJNME]

% % clc
% clearvars
% close all


postTraitement = false;

test = true;
% test = false;

if postTraitement
    setProblem = false;
    solve = false;
    plotResults = true;
    saveParaview = false;
else
    setProblem = true;
    solve = true;
    plotResults = true;
    saveParaview = false;
end

display = false;

%% Model

% solver = 'HistoryField'; % 'HistoryField' or 'BoundConstrainedOptim'
% PFmodel = 'AnisotropicMiehe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
% alpha = 'AT1'; %'AT1', 'AT2'

filename = append('BenchmarkCompressionTest_',PFmodel,'_', alpha,'_',solver);

if test
    filename = append(filename,'_Test');
    pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest','Test',filename);
else
    pathname = fullfile(GetCodesMatthieuPath,'results','BenchmarkCompressionTest',filename);
end

if ~exist(pathname,'dir')
    mkdir(pathname);
end

if setProblem

    %% Datas 
    % [Nguyen, Yvonnet, Waldmann, He, 2020,IJNME]
    
    DIM = 2;
    
    height = 30e-3; %[m]
    width = 15e-3;
    ep = 1;   
    
    holeHeight = height/2; %hauteur trou
    holeDiameter = 12e-3/2;   %diam trou   
    
    
    %% Création du modèle et maillage  

    if test
        cl = 0.25e-3;
        clH = 0.12e-3;
    else
        cl = l_0/2;
        clH = cl;
    end    
    
    domain = DOMAIN(2,[0,0],[width,height]);
    circle = CIRCLE(width/2,height-holeHeight,holeDiameter/2);
    
    model = gmshdomainwithhole(domain,circle,cl,clH,fullfile(pathname,'gmshDomainWithHoleBenchmark'));
        
    %% Modele a gradient d'endommagement
    
    gc = 1.4;   % Taux de libération d'énergie critique (ou ténacité à la rupture) [N/m]
    l_0 = 0.12e-3; % Paramètre de régularisation (largeur de la fissure)
    
    g = @(d) (1-d).^2;  % Fonction degradation energetique    
    
    PF_problem = PF_Problem(solver, PFmodel, alpha, gc, l_0);
    
    % Materiau
    mat_phase = FOUR_ISOT('k',PF_problem.k,'r',PF_problem.r(0));
    mat_phase = setnumber(mat_phase,1);

    % Attribu le matériau
    S_phase = setmaterial(model,mat_phase);
    
    % ici pas besoin de de renseigner des conditions limite car pas de fissure preexitante
    S_phase = final(S_phase);
    
    PF_problem.S_phase = S_phase;
    
    d = calc_init_dirichlet(S_phase);
    
    
    %% Modele deformation elastique
    
    RHO = 1;    % Densité kg/m3
    k = 1e-10;  % Petite rigidité résiduelle artificielle
    E = 12e9;
    NU = 0.3;
            
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
    
    mat = setnumber(mat,1);
    
    S = setoption(model,'DEFO');
    S = setmaterial(S,mat);
    
    PF_problem.S = S;
    
    %% Conditions en déplacement
    
    BLower = LIGNE([0, 0],[width, 0]);
    BUpper = LIGNE([0, height],[width, height]);
    Bgauche = LIGNE([0, 0], [0, height]); 
    Bdroite = LIGNE([width, 0], [width, height]);   
    
    [~, loadNodes] = intersect(S,BUpper);
    
    % Get the nodes on the egdes
    [~,noeudsBLower] = intersect(S, BLower);
    [~,noeudsBUpper] = intersect(S, BUpper);
    [~,noeudsBGauche] = intersect(S, Bgauche);
    [~,noeudsBDroite] = intersect(S, Bdroite);
    noeudsDuBord = [noeudsBGauche; noeudsBDroite; noeudsBLower; noeudsBUpper];
    
    S = final(S);
    PF_problem.DirichletBoundaryConditions{1} = {BLower, 'UY', 0};
    PF_problem.DirichletBoundaryConditions{2} = {[0,0], 'UX', 0};
    
    S = ApplyDirichletBoundaryConditions(S, PF_problem.DirichletBoundaryConditions);
    
    % plot(S,'numnode')
    
    %% Chargement
    
    uIncBig = 8e-8;   %[m] step d<=0.6
    uIncSmall = 2e-8;   % step d>0.6
    maxDep = 25e-6; % 25 µm

    if test
        % maxDep = 14e-6;        
        uIncBig = 16e-8;
        uIncSmall = 4e-8;
    end
    
    save(fullfile(pathname,'problem.mat'), 'PF_problem');
else
    load(fullfile(pathname,'problem.mat'));    
    S = PF_problem.S;
    S_phase = PF_problem.S_phase;
end

%% Resolution
if solve

    d = calc_init_dirichlet(S_phase);
    u = calc_init_dirichlet(S);
    H = calc_energyint(S,u,'positive','intorder','mass');
    
    if display
        PF_problem
        fprintf('\n')
        fprintf('+--------+--------+----------+----------+---------------+\n');
        fprintf('|  Iter  | u [µm] |  min(d)  |  max(d)  |   t [h:m:s]   |\n');
        fprintf('+--------+--------+----------+----------+---------------+\n');
    end
    
    ud=0;
    
    ccBord = 0; % Counter to count the number of times the edge is damaged
    iterBord = 50; % Number of iterations allowed after edge damage 
        
    i=0;    
    while ud < maxDep
        
        i = i+1;
        
        if any(d > 0.6)
            uInc = uIncSmall;
        else
            uInc = uIncBig;
        end        
        
        % inc load
        ud = ud + uInc;
        ud_t(i) = -ud;
        PF_problem.DirichletBoundaryConditions{3} = {loadNodes, 'UY', -ud};
        
        % resolution
        tSolve = tic;
        [u,d,A,H] =  PhaseFieldSolver(S_phase, S, i, PF_problem, H, u, d);
        
        % Force field
        numddl = findddl(S,'UY',loadNodes);
        f = -A(numddl,:)*u;
        f = sum(f);
        
        % calculation of remaining time
        resolutionTime(i) = toc(tSolve);
        
        iterRestant = abs(maxDep-ud)/uInc;
        tempsRestant = resolutionTime(i)*iterRestant;
        
        [hour, minutes, segondes] = GetTime(tempsRestant);        
        
        % Update fields
        dt{i} = d;
        ut{i} = u;
        ft(i) = f;
        Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);

        % Check if the edge is damaged  
        maxdBord = round(max(d(noeudsDuBord)),3);
        if maxdBord >= 1
            ccBord = ccBord + 1;
            fprintf('\n Damaged edge : %d / %d  \n', ccBord, iterBord)
            if ccBord == iterBord
                break
            end
        end
    
        if display
            fprintf('|  %4d  |  %4.2f  | %4f | %4f | %2dh:%2.fm:%2.fs   | %4f \n', ...
                i, ud*1e6, abs(min(dt{i})), abs(max(dt{i})), hour, minutes, segondes, resolutionTime(i));
        end
        
    end

    if display
        fprintf('+--------+--------+----------+----------+---------------+\n');
    end

    resume = PF_problem.resume(sum(resolutionTime));
    
    fprintf(resume)
    
    save(fullfile(pathname,'solution.mat'),'ud_t','dt','ut','Ht','ft','resolutionTime');

else
    load(fullfile(pathname,'solution.mat'));
    resume = PF_problem.resume(sum(resolutionTime));
    fprintf(resume)
end



if plotResults
    
    figure
    plot(ud_t*1e6,-ft*1e-6,'LineWidth',1)
    grid on
    xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
    ylabel("Load in kN/mm",'interpreter','Latex','fontsize',15)
    
    saveas(gcf, fullfile(pathname, 'displacement.png'))

    plotSolution(S_phase,dt{length(dt)});
    saveas(gcf, fullfile(pathname, 'damage.png'))

    figure
    plot(-ud_t*1e6,resolutionTime,'LineWidth',1)
    grid on
    xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
    ylabel("Resolution time in s",'interpreter','Latex','fontsize',15)    
    saveas(gcf, fullfile(pathname, 'time.png'))
    
    % norm dInc
    normDt(1)=0;    
    for i=1:length(dt)-1        
        normDt(i+1) = sqrt(sum((dt{i+1}-dt{i}).^2));
    end    
    figure
    semilogy(-ud_t*1e6, normDt,'LineWidth',1)
    xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
    ylabel("$\Vert d_{i+1} - d_{i} \Vert$",'interpreter','Latex','fontsize',15)
    grid on
    saveas(gcf, fullfile(pathname, 'convergence.png'))

    fprintf('\n plot saved \n')

end

%% Save solutions

calcStress = false;

if saveParaview
    fprintf('\n Paraview \n',i)
    for i=1:length(dt)
        
        text = fprintf('%i / %i \n',i,length(dt));

        di = dt{i};
        ui = ut{i};
        
        if calcStress % Calcul des contraintes
            % update mat
            mats = MATERIALS(S);
            for m=1:length(mats)
                mats{m} = setparam(mats{m},'d',di);
                mats{m} = setparam(mats{m},'u',ui);
            end
            S = actualisematerials(S,mats);
            
            % Stress
            stress = double(mean(calc_sigma(S,ui),4));
            
            Sxx = reshape(stress(1,1,:),[getnbelem(S),1]);
            Syy = reshape(stress(2,1,:),[getnbelem(S),1]);
            Sxy = reshape(stress(3,1,:),[getnbelem(S),1]);       
    
            Svm = sqrt(Sxx.^2+Syy.^2-(Sxx.*Syy)+3*Sxy.^2);
            
            write_vtk_mesh(S,{di,ui},{Sxx,Syy,Sxy,Svm},...
                {'damage','displacement'},{'Sxx','Syy','Sxy','Svm'},...
                pathname,'solution',1,i-1);
        else
            write_vtk_mesh(S,{di,ui},{},{'damage','displacement'},{},pathname,'solution',1,i-1);
        end
        
        if not(i == length(dt))
            fprintf(repmat('\b',1,text));
        end        
        
    end
    make_pvd_file(pathname,'solution',1,length(dt));
end
