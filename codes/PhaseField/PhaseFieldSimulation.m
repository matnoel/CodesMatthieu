function phaseFieldSolution = PhaseFieldSimulation(PFM, displacement,uniteAffichage,signe, display,noeudsDuBord)
% Simulation for a phase field problem
% PFM : phase field model
% displacement : must be a Timemodel or a struct containing :
%     inc0 = displacement.inc0;
%     inc1 = displacement.inc1;
%     umax = displacement.umax;
%     dthreshold = displacement.dthreshold;

%% Load variables
S = PFM.S;
S_phase = PFM.S_phase;

if signe ~= 1 && signe ~= -1
    error("Signe doit etre -1 ou 1")
end

switch uniteAffichage
    case "m"
        coef = 1;
    case "mm"
        coef = 1e3;
    case "µm"
        coef = 1e6;
    otherwise
        error("unité pas exsistente")
end

if isstruct(displacement)
    inc0 = displacement.inc0;
    inc1 = displacement.inc1;
    umax = displacement.umax;
    dthreshold = displacement.dthreshold;

    thresholdSimulation = true;
else
    thresholdSimulation = false;

%     inc0 = linspace(dt0,nt0*dt0,nt0);
%     inc1 = linspace(t0(end)+dt1,t0(end)+nt1*dt1,nt1);
%     t = [inc0,inc1];     
%     T = TIMEMODEL(t);

    displacement = gett(displacement);
    umax = displacement(length(displacement));
end

%% Init

d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
H = calc_energyint(S,u,'positive','intorder','mass');

if display
    fprintf(PFM.resume)
    fprintf('\n')
    fprintf('+--------+--------+----------+----------+-------------+\n');
    fprintf('|  Iter  | u [%s] |  min(d)  |  max(d)  |  t [h:m:s]  |\n',uniteAffichage);
    fprintf('+--------+--------+----------+----------+-------------+\n');
end

ccBord = 0; % Counter to count the number of times the edge is damaged
iterBord = 50; % Number of iterations allowed after edge damage 

ud=0;

i=0;

while ud < umax
    
    i=i+1;
    
    %update displacement conditions
    if thresholdSimulation        
        if any(d > dthreshold)
            uInc = inc1;
        else
            uInc = inc0;
        end        
        ud = ud + uInc;
        iterRestant = abs(umax-ud)/uInc;
    else
        ud = displacement(i);
        % calc uInc for estimation time
        if i>1
            uinc = ud-displacement(i-1);
        else
            uInc = 0;
        end
        iterRestant = length(displacement)-i;
        
    end
    udt(i) = signe*ud;

    %update the dirchelet boundary conditions
    nbCond = length(PFM.DirichletBoundaryConditions); % get the number of boundary conditions
    loadNodes = PFM.DirichletBoundaryConditions{nbCond}{1};
    directions = PFM.DirichletBoundaryConditions{nbCond}{2};
    PFM.DirichletBoundaryConditions{nbCond} = {loadNodes, directions, udt(i)}; %update
    
    % resolution
    tSolve = tic;
    % [u,d,A,H,PFM] =  PhaseFieldSolver(PFM, H, u, d);
    [u,d,A,H,PFM,iter] =  PhaseFieldSolver(PFM, H, u, d);
    
    % Get models back
    S = PFM.S;
    S_phase = PFM.S_phase;

    % Force field
    numddl = findddl(S,directions,loadNodes);
    f = -A(numddl,:)*u;
    f = sum(f);
    
    % calculation of remaining time
    resolutionTime(i) = toc(tSolve);
    
    tempsRestant = resolutionTime(i)*iterRestant;
    
    temps = GetTime(tempsRestant);        
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);

    % Check if the edge is damaged
    maxdBord = round(max(d(noeudsDuBord)),2);
    if maxdBord >= 1
        ccBord = ccBord + 1;
        if display
            fprintf('\n Damaged edge : %d / %d  \n', ccBord, iterBord)
        end            
        if ccBord == iterBord
            break
        end
    end

    if display
        dep = ud*coef;
        maxd = abs(max(dt{i}));
        mind = abs(min(dt{i}));

        fprintf('|  %4d  |  %4.2f  | %4.3f | %4.3f | %s  | %4.3f | %d \n', ...
            i, dep, mind, maxd, temps, resolutionTime(i), iter);
    end
    
    if iter>60
        fprintf("Arret car on arrive pas a converger iter = %s \n",iter)
        break
    end

end

if display
    fprintf('+--------+--------+----------+----------+-------------+\n');
end

fprintf("\n"+GetTime(sum(resolutionTime))+"")

phaseFieldSolution = PhaseFieldSolution(PFM,udt,dt,ut,Ht,ft,resolutionTime);

