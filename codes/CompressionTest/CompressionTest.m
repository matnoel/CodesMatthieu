%% Compression test FCBA %%
clc;    clearvars;  close all;

solve = false;

plot = true;
saveParaview = false;

optimesh = true;
test = false;

materiau='ELAS_ISOT_TRANS'; % 'ELAS_ISOT', 'ELAS_ISOT_TRANS' 
option = 'CONT'; % 'DEFO' plane strain, 'CONT' plane stress  
PFmodel = 'AnisotropicHe'; % 'Isotropic', 'AnisotropicAmor', 'AnisotropicMiehe', 'AnisotropicHe'
alpha = 'AT1'; %'AT1', 'AT2'

% Builds the name and location of the file
filename = append(materiau,'_',option,'_',PFmodel,'_', alpha);

if optimesh
    filename = append(filename,'_OptiMesh');
end

if test
    filename = append(filename,'_Test');
    pathname = fullfile(GetCodesMatthieuPath,'results','CompressionTest','Test',filename);
else
    pathname = fullfile(GetCodesMatthieuPath,'results','CompressionTest',filename);
end

% Create the file if it does not exist
if ~exist(pathname,'dir')
    mkdir(pathname);
end

display = true;

%% Datas

DIM = 2;

height = 12e-2; %[m]
width = 9e-2;
ep = 2e-2;   

holeHeight = 3.5e-2; 
holeDiameter = 1e-2;

coefZone = 3;

%% Model creation and meshing

l_0 = width/100; % Regulation parameter (crack width)

if test
    l_0 = width/50;   
end

domain = DOMAIN(DIM,[-width/2,0],[width/2,height]);
circle = CIRCLE(0.0,height-holeHeight,holeDiameter/2);

model = ConstruitModelEssaiCompression(domain, circle, l_0, coefZone, PFmodel, optimesh, pathname);

%% Damage gradient model

gc = 1.4;   % Critical energy release rate (or fracture toughness) 
k = 1e-10;  % Small artificial residual stiffness 
H = 0;  % Internal energy 

% Material
if alpha == 'AT1'
    mat_phase = FOUR_ISOT('k',3/4*gc*l_0,'r',2*H);  % AT1
    constAT1 = 3*gc/8/l_0;
else
    mat_phase = FOUR_ISOT('k',gc*l_0,'r',gc/l_0+2*H);   % AT2
    constAT1 = 0;
end
mat_phase = setnumber(mat_phase,1);

% Allocate material
S_phase = setmaterial(model,mat_phase);

% here no need to enter boundary conditions as no pre-existing cracks
S_phase = final(S_phase);

%% Elastic deformation model

d = calc_init_dirichlet(S_phase);
g = @(d) (1-d).^2;  % Energy degradation function
RHO = 1;

if materiau == "ELAS_ISOT"
    E = 11580e6; % 11580 MPa
    NU = 0.3;
    mat = ELAS_ISOT('E',E,'NU',NU,'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
else
    % Transverse Young modulus
    ET = 500e6; % [Pa]
    % Longitudinal shear modulus
    GL = 175e6; % [Pa]
    % Longitudinal Young modulus
    EL = 11580e6; % [Pa]
    % Longitudinal Poisson ratio
    NUL = 0.44;
    % Transverse Poisson ratio
    NUT = 0.02;
    % Material
    mat = ELAS_ISOT_TRANS('AXISL',[0;1;0],'AXIST',[1;0;0], ...
        'EL',EL,'ET',ET,'NUL',NUL,'NUT',NUT,'GL',GL, ...
        'RHO',RHO,'DIM3',ep,'d',d,'g',g,'k',k,'u',0,'PFM',PFmodel);
end

mat = setnumber(mat,1);
S = setoption(model,option);
S = setmaterial(S,mat);

%% Retrieves nodes to apply boundary conditions

% Geometric object
PointGauche = POINT([-width/2, 0]);
BLower = LIGNE([-width/2, 0],[width/2, 0]);
BUpper = LIGNE([-width/2, height],[width/2, height]);
Bgauche = LIGNE([-width/2, 0], [-width/2, height]); 
Bdroite = LIGNE([width/2, 0], [width/2, height]);   

S = final(S); % Finalize the mesh

S = addcl(S,BLower,{'UX','UY'});
% S = addcl(S,BLower,{'UY'});
% S = addcl(S,PointGauche,{'UX'});


% Retrieves nodes from the edge
[~,noeudsBLower] = intersect(S, BLower);
[~,noeudsBUpper] = intersect(S, BUpper);
[~,noeudsBGauche] = intersect(S, Bgauche);
[~,noeudsBDroite] = intersect(S, Bdroite);

noeudsDuBord = [noeudsBGauche; noeudsBDroite; noeudsBLower; noeudsBUpper];

% Nodes of the circle
[~,nodesBCircle] = intersect(S,circle);  
coord = getcoord(S.node(nodesBCircle));

% nodes of the lower circle
% loadNodes = nodesBCircle(find(coord(:,2)-(height-holeHeight)<1e-6));    
[~, loadNodes] = intersect(S,BUpper);

%% Displacement parameter

if test 
    uInc = 8e-6; % Increment when d<=0.2
    uIncEndo = 2e-6;
else
    uInc = 8e-7; % Increment when d<=0.2
    uIncEndo = 2e-7;
end

ud = 0; % Initial displacement

%% Resolution

% initializes variables
d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
H = calc_energyint(S,u,'positive','intorder','mass');

ccBord = 0; % Counter to count the number of times the edge is damaged
iterBord = 20; % Number of iterations allowed after edge damage

if solve
    
    if display
        fprintf('\n')
        fprintf('+--------+-----------+----------+----------+------------+\n');
        fprintf('|  Iter  |  u [mm]   |  min(d)  |  max(d)  |  norm(dInc) |\n');
        fprintf('+--------+-----------+----------+----------+------------+\n');
    end
    
    i = 1;
    while ccBord < iterBord
        
        % imposed displacement
        if max(d) <= 0.2
            ud = ud - uInc;
        else
            ud = ud - uIncEndo;
        end
        
        d_old = d;

        [u,d,f,H] =  ResolutionCompressionTest(S_phase,S,BLower,loadNodes, H, u, ud, constAT1);        
        
        % Update fields -------------------------------------------------------
        dt{i} = d;
        ut{i} = u;
        udt(i) = ud;
        ft(i) = f;
        Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);
        
        % Convergence verification
        dInc = d - d_old;
        norm = sqrt(sum(dInc.^2));
    
        if display
            fprintf('| %4d   | %6.3e | %8.2e | %8.2e |  %8.3e  |\n',i, abs(ud)*1e3, min(abs(dt{i})), max(abs(dt{i})), norm);
        end
        
        % Check if the edge is damaged  
        maxdBord = max(round(d(noeudsDuBord),2));        
        if maxdBord >= 1
            ccBord = ccBord + 1;
            fprintf('\n Damaged edge \n')
        end
        i=i+1;    
    end

    if display
        fprintf('+-----------+-----------+----------+----------+------------+\n');
    end

    save(fullfile(pathname,'solution.mat'),'udt','dt','ut','Ht','ft','S','S_phase');
else
    load(fullfile(pathname,'solution.mat'));
end


%% Post treatment 


if plot
    % plotDomain(model,'legend',false);
    % plotModel(model,'Color','k','FaceColor','c','FaceAlpha',1,'legend',false);    
    
    clear plot
    figure
    plot(-udt*1e3,ft*1e-3)
    grid on
    xlabel("Displacement in mm",'interpreter','Latex','fontsize',15)
    ylabel("Load in kN",'interpreter','Latex','fontsize',15)
    
    saveas(gcf, fullfile(pathname, 'displacement.png'))

    plotSolution(S_phase,dt{length(dt)});
    saveas(gcf, fullfile(pathname, 'damage.png'))
    
    
    % norm dInc
    normDt(1)=0;
    for i=1:length(dt)-1
        normDt(i+1) = sqrt(sum((dt{i+1}-dt{i}).^2));
    end
    figure
    semilogy(-udt*1e3, normDt)
    xlabel("Displacement in mm",'interpreter','Latex','fontsize',15)
    ylabel("$\Vert d_{i+1} - d_{i} \Vert$",'interpreter','Latex','fontsize',15)
    grid on
    saveas(gcf, fullfile(pathname, 'convergence.png'))

    fprintf('\n plot saved')
end

%% Save solutions

if saveParaview
    fprintf('\n Paraview \n',i)
    for i=1:length(dt)

        fprintf(' %i / %i \n',i,length(dt))

        di = dt{i};
        ui = ut{i};
        Hi = Ht{i};
           
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
        
        write_vtk_mesh(S,{di,ui},{Hi,Sxx,Syy,Sxy,Svm},...
            {'damage','displacement'},{'internal energy density history','Sxx','Syy','Sxy','Svm'},...
            pathname,'solution',1,i-1);
    end
    make_pvd_file(pathname,'solution',1,length(dt));
    fprintf('\n Paraview saved')
end