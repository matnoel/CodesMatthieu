classdef PhaseFieldSolution
    % [udt,dt,ut,Ht,ft,resolutionTime]
    properties
        PFM % Phase field model
        udt % Imposed displacement
        dt % Damaged fields
        ut % Displacement field
        Ht % Internal energy
        ft % Force related to the imposed displacement
        resolutionTime % Resolution time for each solver iteration
    end

    methods

        function obj = PhaseFieldSolution(PFM,udt,dt,ut,Ht,ft,resolutionTime)
            obj.PFM = PFM;
            obj.udt = udt;
            obj.dt = dt;
            obj.ut = ut;
            obj.Ht = Ht;
            obj.ft = ft;
            obj.resolutionTime = resolutionTime;
        end

        function [] = PlotResults(obj, pathname, unite)

            switch unite
                case "m"
                    coef = 1;
                case "mm"
                    coef = 1e3;
                case "µm"
                    coef = 1e6;
                otherwise
                    error("unité pas exsistente")
            end

            mats = MATERIALS(obj.PFM.S);
            mat = mats{1};
            ep = getparam(mat,'DIM3');
            
            dept = abs(obj.udt*coef);
            forcet = abs(obj.ft/1e6);

            figure
            plot(dept,forcet,'LineWidth',1)
            grid on
            if unite == "µm"
                ytitle = "Displacement in $\mu$m";
%                 ytitle = "D\'eplacement en $\mu$m";
            else
                ytitle = "Displacement in $"+unite+"$";
            end
            
            xlabel(ytitle,'interpreter','Latex','fontsize',15)
            ylabel("Load in kN/mm",'interpreter','Latex','fontsize',15)
%             ylabel("Force de compression en kN/mm",'interpreter','Latex','fontsize',15)
            
            saveas(gcf, fullfile(pathname, 'displacement.png'))
            
            plotSolution(obj.PFM.S_phase,obj.dt{length(obj.dt)});
            saveas(gcf, fullfile(pathname, 'damage.png'))
            
%             figure
%             plot(-obj.udt,obj.resolutionTime,'LineWidth',1)
%             grid on
%             xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
%             ylabel("Resolution time in s",'interpreter','Latex','fontsize',15)    
%             saveas(gcf, fullfile(pathname, 'time.png'))
            
%             % norm dInc
%             normDt(1)=0;    
%             for i=1:length(obj.dt)-1        
%                 normDt(i+1) = sqrt(sum((obj.dt{i+1}-obj.dt{i}).^2));
%             end
%             figure
%             semilogy(-obj.udt*1e6, normDt,'LineWidth',1)
%             xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
%             ylabel("$\Vert d_{i+1} - d_{i} \Vert$",'interpreter','Latex','fontsize',15)
%             grid on
%             saveas(gcf, fullfile(pathname, 'convergence.png'))
            
            fprintf("\n saved \n")
        end

        function [] = SaveParaview(obj, pathname)
            
            S = obj.PFM.S;

            calcStress = false;
            
            fprintf('\n Paraview \n')
            for i=1:length(obj.dt)
                
                text = fprintf('%i / %i \n',i,length(obj.dt));
            
                ui = obj.ut{i};
                di = obj.dt{i};
                
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
                
                if not(i == length(obj.dt))
                    fprintf(repmat('\b',1,text));
                end        
                
            end
            make_pvd_file(pathname,'solution',1,length(obj.dt));
        end

        function [] = MakeMovie(obj, pathname)
            
            options = {'plotiter',true,'plottime',false};
            framerate = 80;
            
            T = TIMEMODEL(obj.udt);

            dt = TIMEMATRIX(obj.dt,T,size(obj.dt));

            evolSolution(obj.PFM.S_phase, dt,'FrameRate',framerate,'filename','damage','pathname',pathname,options{:});

        end
        
    end
end