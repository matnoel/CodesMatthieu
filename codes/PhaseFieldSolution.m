classdef PhaseFieldSolution
    % [udt,dt,ut,Ht,ft,resolutionTime]
    properties
        PFM % Phase field model
        udt % Imposed displacement
        dt % Damaged field
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

        function [] = PlotResults(obj, pathname)
            
            fprintf(obj.PFM.resume)

            figure
            plot(obj.udt*1e6,-obj.ft*1e-6,'LineWidth',1)
            grid on
            xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
            ylabel("Load in kN/mm",'interpreter','Latex','fontsize',15)
            
            saveas(gcf, fullfile(pathname, 'displacement.png'))
            
            plotSolution(obj.PFM.S_phase,obj.dt{length(obj.dt)});
            saveas(gcf, fullfile(pathname, 'damage.png'))
            
            figure
            plot(-obj.udt*1e6,obj.resolutionTime,'LineWidth',1)
            grid on
            xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
            ylabel("Resolution time in s",'interpreter','Latex','fontsize',15)    
            saveas(gcf, fullfile(pathname, 'time.png'))
            
            % % norm dInc
            % normDt(1)=0;    
            % for i=1:length(dt)-1        
            %     normDt(i+1) = sqrt(sum((dt{i+1}-dt{i}).^2));
            % end    
            % figure
            % semilogy(-ud_t*1e6, normDt,'LineWidth',1)
            % xlabel("Displacement in $\mu m$",'interpreter','Latex','fontsize',15)
            % ylabel("$\Vert d_{i+1} - d_{i} \Vert$",'interpreter','Latex','fontsize',15)
            % grid on
            % saveas(gcf, fullfile(pathname, 'convergence.png'))
            
            fprintf("\n saved \n")
        end
        
    end
end