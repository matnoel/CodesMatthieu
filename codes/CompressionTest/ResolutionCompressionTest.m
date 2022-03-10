function [u,d,f,H] =  ResolutionCompressionTest(S_phase,S,BLower,loadNodes, H, u, ud, constAT1)

% ------------------------Internal energy field ---------------------------
h_old = getvalue(H);
H = calc_energyint(S,u,'positive','intorder','mass');
h = getvalue(H);
for p=1:getnbgroupelem(S)
    he = double(h{p});
    he_old = double(h_old{p});
    rep = find(he <= he_old);
    he(rep) = he_old(rep);
    h{p} = he;
end
H = FEELEMFIELD(h,'storage',getstorage(H),'type',gettype(H),'ddl',getddl(H));

% ------------------------------Phase field -------------------------------
mats_phase = MATERIALS(S_phase);
for m=1:length(mats_phase)
    r = getparam(mats_phase{m},'r');

    if constAT1 == 0
        r = r + 2*H;        
    else
        r = 2*H;
    end

    mats_phase{m} = setparam(mats_phase{m},'r',r);
    
end


S_phase = actualisematerials(S_phase,mats_phase);

[A_phase,b_phase] = calc_rigi(S_phase);

if constAT1 == 0
    b_phase = -b_phase + bodyload(S_phase,[],'QN',2*H); %AT2
else
    f = 2*H-constAT1;
    fPlus = (abs(f)+f)/2;
    b_phase = -b_phase + bodyload(S_phase,[],'QN',fPlus); %AT1
end

d = A_phase\b_phase;
d = unfreevector(S_phase,d);

% --------------------------Displacement field ----------------------------
mats = MATERIALS(S);
for m=1:length(mats)
    mats{m} = setparam(mats{m},'d',d);
    mats{m} = setparam(mats{m},'u',u);
end
S = actualisematerials(S,mats);
S = removebc(S);

S = addcl(S,loadNodes,'UY',ud);

S = addcl(S,BLower,{'UX','UY'});
% S = addcl(S,BLower,{'UY'});
% S = addcl(S,1,{'UX'});


[A,b] = calc_rigi(S,'nofree');
b = -b;

u = freematrix(S,A)\b;
u = unfreevector(S,u);

% -----------------------------Force Field --------------------------------
numddl = findddl(S,'UY',loadNodes);
f = -A(numddl,:)*u;
f = sum(f);


