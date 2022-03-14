function phaseFieldSolution = Treshold_PhaseFieldSimulation(PFM, loadNodes,smallInc, bigInc, maxDep, treshold, display)

%% Load variables
S = PFM.S;
S_phase = PFM.S_phase;

%% Init

d = calc_init_dirichlet(S_phase);
u = calc_init_dirichlet(S);
H = calc_energyint(S,u,'positive','intorder','mass');

if display
    PFM
    fprintf('\n')
    fprintf('+--------+--------+----------+----------+-------------+\n');
    fprintf('|  Iter  | u [µm] |  min(d)  |  max(d)  |  t [h:m:s]  |\n');
    fprintf('+--------+--------+----------+----------+-------------+\n');
end

ud=0;

ccBord = 0; % Counter to count the number of times the edge is damaged
iterBord = 50; % Number of iterations allowed after edge damage 
    
i=0;    
while ud < maxDep
    
    i = i+1;
    
    if any(d > treshold)
        uInc = smallInc;
    else
        uInc = bigInc;
    end        
    
    % inc load
    ud = ud + uInc;
    udt(i) = -ud;
    PFM.DirichletBoundaryConditions{3} = {loadNodes, 'UY', -ud};
    
    % resolution
    tSolve = tic;
    [u,d,A,H,PFM] =  PhaseFieldSolver(PFM, H, u, d);
    
    % Get models back
    S = PFM.S;
    S_phase = PFM.S_phase;

    % Force field
    numddl = findddl(S,'UY',loadNodes);
    f = -A(numddl,:)*u;
    f = sum(f);
    
    % calculation of remaining time
    resolutionTime(i) = toc(tSolve);
    
    iterRestant = abs(maxDep-ud)/uInc;
    tempsRestant = resolutionTime(i)*iterRestant;
    
    temps = GetTime(tempsRestant);        
    
    % Update fields
    dt{i} = d;
    ut{i} = u;
    ft(i) = f;
    Ht{i} = reshape(double(mean(H,4)),[getnbelem(S),1]);

%         % Check if the edge is damaged  
%         maxdBord = round(max(d(noeudsDuBord)),3);
%         if maxdBord >= 1
%             ccBord = ccBord + 1;
%             if display
%                 fprintf('\n Damaged edge : %d / %d  \n', ccBord, iterBord)
%             end            
%             if ccBord == iterBord
%                 break
%             end
%         end

    if display
        fprintf('|  %4d  |  %4.2f  | %4f | %4f | %s  | %4f \n', ...
            i, ud*1e6, abs(min(dt{i})), abs(max(dt{i})), temps, resolutionTime(i));
    end
    
end

resume = PFM.resume;

if display
    fprintf('+--------+--------+----------+----------+-------------+\n');
    fprintf(resume)
end

fprintf("\n"+GetTime(sum(resolutionTime))+"")

phaseFieldSolution = PhaseFieldSolution(PFM,udt,dt,ut,Ht,ft,resolutionTime);

