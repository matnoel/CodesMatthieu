function [u,d,A,H,PFM,iter] = PhaseFieldSolver(PFM, H, u, d)

%% ----------------------------Import variables ---------------------------

solver = PFM.solver; % 'HistoryField' or 'BoundConstrainedOptim'

S = PFM.S;
S_phase = PFM.S_phase;

if strcmpi(solver,'BoundConstrainedOptim')
    optimFun = PFM.optimFun;
    options = PFM.options;
end



%% ---------------- Resolution avec gestion de convergence ----------------

% ici on va vérifier quon a converger au sens de bourdin 
% [B. Bourdin, G.A. Francfort, J-J. Marigo, Numerical experiments in revisited brittle fracture, 2000]
%  condition = max(abs(d-dold))<tol
tol=1e-2;
condition=false;
iter=0;

dold=d;

while ~condition
    iter = iter + 1;

    %% ------------------------- Internal energy field ------------------------
    
    switch solver
        case 'HistoryField'
            h_old = getvalue(H);
            H = calc_energyint(S,u,'positive','intorder','mass');
            h = getvalue(H);
            for p=1:getnbgroupelem(S)
                he = double(h{p});
                he_old = double(h_old{p});
                rep = find(he <= he_old);
                he(rep) = he_old(rep);
                h{p} = MYDOUBLEND(he);
            end
            H = FEELEMFIELD(h,'storage',getstorage(H),'type',gettype(H),'ddl',getddl(H));
        otherwise
            H = calc_energyint(S,u,'positive','intorder','mass');
    end
    
    %% ----------------------------- Phase field ------------------------------
    
    mats_phase = MATERIALS(S_phase);
    for m=1:length(mats_phase)
       mats_phase{m} = setparam(mats_phase{m},'r',PFM.r(H{m}));    
    end
    S_phase = actualisematerials(S_phase,mats_phase);
    
    [A_phase,b_phase] = calc_rigi(S_phase);
    
    b_phase = -b_phase + bodyload(S_phase,[],'QN',PFM.F(H));
    
    % Résolution d
    if any(b_phase>0)
        switch solver
            case 'HistoryField'
                d = A_phase\b_phase;            
            otherwise
                d0 = freevector(S_phase,d);
                lb = d0;
                lb(lb==1) = 1-eps;
                ub = ones(size(d0));
                switch optimFun
                    case 'lsqlin'
                        [d,err,~,exitflag,output] = lsqlin(A_phase,b_phase,[],[],[],[],lb,ub,d0,options);
                    case 'lsqnonlin'
                        fun = @(d) funlsqnonlinPF(d,A_phase,b_phase);
                        [d,err,~,exitflag,output] = lsqnonlin(fun,d0,lb,ub,options);
                    case 'fmincon'
                        fun = @(d) funoptimPF(d,A_phase,b_phase);
                        [d,err,exitflag,output] = fmincon(fun,d0+eps,[],[],[],[],lb,ub,[],options);
                end
        end 
        d = unfreevector(S_phase,d);
    end
    
    %% ------------------------ Displacement field ----------------------------
    
    % Update mat
    mats = MATERIALS(S);
    for m=1:length(mats)
        mats{m} = setparam(mats{m},'d',d);
        mats{m} = setparam(mats{m},'u',u);
    end
    S = actualisematerials(S,mats);
    S = removebc(S);
    S = ApplyDirichletBoundaryConditions(S, PFM.DirichletBoundaryConditions);
    
    [A,b] = calc_rigi(S,'nofree');
    b = -b;
    
    u = freematrix(S,A)\b;
    u = unfreevector(S,u);

    dincMax=max(abs(d-dold));
%     max(d)
%     max(dold)
    condition = dincMax < tol;

    dold=d;

    if iter>60
        break
    end

end

%% ---------------------------- Updtate PFM -------------------------------

PFM.S = S;
PFM.S_phase = S_phase;