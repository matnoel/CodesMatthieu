function [u,d,A,H] = PhaseFieldSolver(S_phase, S, iter, PF_Problem, H, u, d)

%% ----------------------------Import variables ---------------------------

solver = PF_Problem.solver; % 'HistoryField' or 'BoundConstrainedOptim'

if strcmpi(solver,'BoundConstrainedOptim')
    optimFun = PF_Problem.optimFun;
    options = PF_Problem.options;
end


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
   mats_phase{m} = setparam(mats_phase{m},'r',PF_Problem.r(H{m}));    
end
S_phase = actualisematerials(S_phase,mats_phase);

[A_phase,b_phase] = calc_rigi(S_phase);
b_phase = -b_phase + bodyload(S_phase,[],'QN',PF_Problem.F(H));

if any(b_phase>0)

% Resolution
% d_old = d;
switch solver
    case 'HistoryField'
        d = A_phase\b_phase;            
    otherwise
        d0 = freevector(S_phase,d);
        ub = ones(size(d0));
        if isequal(d0,ub)
            d = d0;
        else
            lb = d0;
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
end 

d = unfreevector(S_phase,d);
% max(d-d_old)
end



%% ------------------------ Displacement field ----------------------------

mats = MATERIALS(S);
for m=1:length(mats)
    mats{m} = setparam(mats{m},'d',d);
    mats{m} = setparam(mats{m},'u',u);
end
S = actualisematerials(S,mats);
S = removebc(S);
S = ApplyDirichletBoundaryConditions(S, PF_Problem.DirichletBoundaryConditions);

[A,b] = calc_rigi(S,'nofree');
b = -b;

u = freematrix(S,A)\b;
u = unfreevector(S,u);

