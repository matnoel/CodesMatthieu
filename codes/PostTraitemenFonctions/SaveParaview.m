function [] = SaveParaview(ut, dt, PFM)

calcStress = false;

fprintf('\n Paraview \n',i)
for i=1:length(dt)
    
    text = fprintf('%i / %i \n',i,length(dt));

    ui = ut{i};
    di = dt{i};
    
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