function S = ApplyDirichletBoundaryConditions(S, DirichletBoundaryConditions)

if isempty(DirichletBoundaryConditions)
    return
end

S = removebc(S);

for i=1:length(DirichletBoundaryConditions)
    
    nodes = DirichletBoundaryConditions{i}{1};
    directions = DirichletBoundaryConditions{i}{2};
    values = DirichletBoundaryConditions{i}{3};
        
    if values == 0
        S = addcl(S,nodes,directions);
    else
        S = addcl(S,nodes,directions, values);
    end

    

end