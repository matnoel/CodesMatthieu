classdef PF_Problem

    properties
        PFmodel % Isotropic, AnisotropicAmor, AnisotropicMiehe, AnisotropicHe
        model % AT1 or AT2
        gc % Critical energy release rate
        l_0 % Regularization length
        S % Displacement model
        S_phase % Phase field model
        DirichletBoundaryConditions % Dirichlet Boundary Conditions

        solver % HistoryField or BoundConstrainedOptim
        optimFun % 'fmincon' or 'lsqnonlin'
        options % options
    end

    methods

        function obj = PF_Problem(solver, PFmodel, model, gc, l_0)

            obj.solver = solver;
            obj.PFmodel = PFmodel;

            if ~strcmpi(model, 'AT1') && ~strcmpi(model, 'AT2')
                error('Wrong regulization model. Must be AT1 or AT2')
            end

            obj.model = model;
            obj.gc = gc;
            obj.l_0 = l_0;

            if strcmpi(solver,'BoundConstrainedOptim')
                [obj.optimFun, obj.options] = Get_BoundConstrainedOptim_options(obj);
            end            
        end

        function k = k(obj)
            % Diffusion parameter
            switch obj.model
                case 'AT1'
                    k = 3/4*obj.gc*obj.l_0;
                case 'AT2'
                    k = obj.gc*obj.l_0;
            end
        end

        function r = r(obj, H)
            % Reaction parameter
            switch obj.model
                case 'AT1'
                    r = 2*H;
                case 'AT2'
                    r = 2*H + obj.gc/obj.l_0;
            end
        end

        function F = F(obj, H)
            % Source therme 
            switch obj.model
                case 'AT1'
                    f = 2*H - (3*obj.gc/8/obj.l_0);
                    F = (f+abs(f))/2;    
                case 'AT2'
                    F = 2*H;
            end
        end

        function resume = resume(obj, time)
            
            DIM = getdim(getgroupelem(obj.S,1));

            [hour, minutes, segondes] = GetTime(time);

            resume = "\n dim : "+DIM+"" + ...
                    "\n solver : "+obj.solver+"" + ...
                    "\n PFmodel : "+obj.PFmodel+"" + ...
                    "\n model : "+obj.model+"" + ...
                    "\n gc : "+obj.gc+"" + ...
                    "\n l_0 : "+obj.l_0+"" + ...
                    "\n nb elements : "+getnbelem(obj.S)+"" + ...
                    "\n nb nodes : "+getnbnode(obj.S)+"" + ...
                    "\n nb dofs : "+getnbnode(obj.S)*DIM+"" + ...
                    "\n elapsed time : "+hour+"h:"+minutes+"m:"+segondes+"s \n";
        end        
      
    end

    methods (Access = private)
        function [optimFun, options] = Get_BoundConstrainedOptim_options(obj)
            if strcmpi(obj.solver,'BoundConstrainedOptim')
                optimFun = 'lsqlin'; % 'lsqlin', 'fmincon' or 'lsqnonlin'
                % optimFun = 'lsqnonlin'; 
                % optimFun = 'fmincon';
            
                displayoptim = 'off';
                % displayoptim = 'iter';
                % displayoptim = 'iter-detailed';
                % displayoptim = 'final';
                % displayoptim = 'final-detailed';
            
                % tolX = 1e-6; % tolerance on the parameter value 1e-6
                % tolFun = 1e-6; % tolerance on the function value 1e-6
                % maxFunEvals = Inf; % maximum number of function evaluations
            
                % optimAlgo = 'interior-point';
                optimAlgo = 'trust-region-reflective';
                % optimAlgo = 'sqp';
                % optimAlgo = 'active-set';
                % optimAlgo = 'levenberg-marquardt';
            
                % options  = optimoptions(optimFun,'Display',displayoptim,'StepTolerance',tolX,'FunctionTolerance',tolFun,...
                %     'OptimalityTolerance',tolFun...%,'MaxFunctionEvaluations',maxFunEvals...%,'Algorithm',optimAlgo...
                %     ,'SpecifyObjectiveGradient',true...
                %     );
                % options  = optimoptions(optimFun,'Display',displayoptim,'SpecifyObjectiveGradient',true);
                
                % options  = optimoptions(optimFun,'Display',displayoptim, ...
                %    'StepTolerance', tolX,'FunctionTolerance',tolFun,'OptimalityTolerance',tolFun,...
                %    'SpecifyObjectiveGradient',true);

                options = optimoptions(optimFun,'Display',displayoptim,'Algorithm',optimAlgo);
            end
        end
    end
end