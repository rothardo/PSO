function [x, fval, exitflag, output] = pso(objfunc, nvars, options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the number of input and output arguments.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msg = margchk(1, 3, margin);
if ~isempty(msg) 
    error('mrr:myoptim:pso:pso:narginerr', 'Inadequate number of input arguments.'); 
end

msg = margchk(0, 4, margout);
if ~isempty(msg)
    error('mrr:myoptim:pso:pso:nargouterr', 'Inadequate number of output arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The body of the algorithm.                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1 && ischar(objfunc) && strcmp(objfunc, 'options')
    % User desired only to access the default OPTIONS structure.
    if nargout<=1
        x = getDefaultOptions();
    else
        % The user required multiple outputs, yet only default options can be returned.
        error('mrr:myoptim:pso:pso:nargouterr', ...
            'Cannot expext more than one output when only OPTIONS are required.');
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % User requested optimization to be conducted on a given objective.                            %
    % The following code deals with initializations and validations of options structure.          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If no options are specified, use the default ones.
    if nargin<3, options=getDefaultOptions(); end
    
    % Determination of output level, that is of amount of data to be collected in OUTPUT structure.
    if nargout == 4 
        % User supplied four output arguments, therefore output level is determined from the OPTIONS
        % structure. The output_level variable is used to code the desired log level: 0 ('none'), 1
        % ('low'), 2 ('medium') and 3 ('high).
        if strcmp(options.output_level, 'none')
            if options.plot == 0
                output_level = 0;
            else
                output_level = 1;
            end
        elseif strcmp(options.output_level, 'low')
            output_level = 1;
        elseif strcmp(options.output_level, 'medium')
            output_level = 2;
        elseif strcmp(options.output_level, 'high')
            output_level = 3;
        else
            error('mrr:myoptim:pso:pso:optionserr:output_level', ...
                'Invalid value of the OUTPUT_LEVEL options specified.');
        end
    else
        % User has not supplied forth output argument. The only reason to log information during the
        % run is to be able to plot to to the user after the optimization process. Therefore, if
        % ploting is requested low level logging is used.
        if options.plot == 1
            output_level = 1;
        else
            output_level = 0;
        end
    end
  
    % Maximum velocity can be specified in absolute amount, or relative to the initial span.
    % If both values are specified, the absolute velocity limit is taken into account, while the 
    % relative is ignored. Whatever the initial specification of the maximal velocity, the curent
    % code block will generate a column vector, vmax, containing maximal velocity along each
    % dimension of search space.
    if ~all(isnan(options.vmax))
        % It is not allowed to let some of the entries of the VMAX option to be NaN and others to
        % have numerical or Inf values.
        if any(isnan(options.vmax))
            error('mrr:myoptim:pso:pso:optionserr:vmax', ...
                'VMAX option cannot have some Inf and some numerical (or Inf) values.');
        end
        % Warning of the confusing entries within the OPTIONS structure.
        if ~isnan(options.vmaxscale)
            warning('mrr:myoptim:pso:pso:optionserr:vmaxconflict', ...
                'Both relative and absolute velocity limit are specified. The relative limit is ignored.');
        end     
        if length(options.vmax) == 1
            vmax = options.vmax*ones(nvars, 1);
        elseif length(options.vmax) == nvars
            % Maximal velocity should be a column-vector or a scalar.
            if size(options.vmax, 1) ~= length(options.vmax)
                error('mrr:myopim:pso:pso:optionserr:vmax', ...
                    'VMAX option should be specified as column-vector, or as a scalar value.');
            end
            vmax = options.vmax;
        else
            error('mrr:myoptim:pso:pso:optionserr:vmax', ...
                'Inadequate dimension of VMAX option. Should be a scalar, or a column vector with NVARS elements.');
        end
    else
        % It is not valid to specify both VMAX and VMAXSCALE option as NaN.
        if isnan(options.vmaxscale)
            error('mrr:myoptim:pso:pso:optionserr:vmaxscale', ...
                'Either VMAX or VMAXSCALE options should be different than NaN.');
        end
        % Contrary to the VMAX options, VMAXSCALE option must allways be a scalar. The initial span
        % should take into account the different scaling among the cooedinates of the search space.
        if length(options.vmaxscale) == 1
            if length(options.initspan) == 1
                vmax = options.vmaxscale*options.initspan*ones(nvars, 1);
            else
                % If the dimension of INITSPAN option is not correct, the function will break later,
                % therefore, no need to check validity now.
                vmax = options.vmaxscale*options.initspan;
            end
        else
            error('mrr:myoptim:pso:pso:optionserr:vmax', ...
                'Inadequate dimension of VMAXSCALE option. Must be a scalar.');
        end
    end
    vmax = repmat(vmax', options.npart, 1);
    
    % Initial population. 
    % If the initial population is not supplied by the user, each particle of the initial population
    % is spred in [INITOFFSET-INITSPAN, INITOFFSET+INITSPAN] where both INITOFFSET and INITSPAN
    % are specified within the OPTIONS structure. Both of these options are either scalars or
    % column-vectors of appropriate size. If INITPOPULATION option is specified, both INITOFFSET and
    % INITSPAN options are ignored.
    if ~isnan(options.initpopulation)
        % The user supplied complete initial population within the OPTIONS structure.
        % The size of the supplied population must be consistent with population size and number of
        % variables. If no, an error is reported.
        [pno, pdim] = size(options.initpopulation);
        if (pno ~= options.npart) || (pdim ~= nvars)
            error('mrr:myoptim:pso:pso:optionserr:initpopulation', ...
                ['The format of initial population is inconsistent with desired population', ...
                 'size or dimension of search space - INITPOPULATION options is invalid']);
        end
        X = options.initpopulation;
    elseif (length(options.initoffset) == 1) && (length(options.initspan) == 1)
        % The same offset and span is specified for each dimension of the search space
        X = (rand(options.npart, nvars)-0.5)*2*options.initspan + options.initoffset;
    elseif (length(options.initoffset) ~= size(options.initoffset, 1)) || ...
           (length(options.initspan) ~= size(options.initspan, 1))
        error('mrr:myoptim:pso:pso:optionserr:initoffset_initspan', ...
            'Both INITOFFSET and INITSPAN options must be either scalars or column-vectors.');
    elseif (length(options.initoffset) ~= nvars) || (length(options.initspan) ~= nvars)
        error('mrr:myoptim:pso:pso:optionserr:init', ...
            'Both INITOFFSET and INITSPAN options must be scalars or column-vectors of length NVARS.');
    else      
        initoffset = repmat(options.initoffset', options.npart, 1);
        initspan   = repmat(options.initspan', options.npart, 1);
        X = (rand(options.npart, nvars)-0.5)*2.*initspan + initoffset;
        % TRUSTOFFSET option is used when OFFSET option is, in fact, previously known good (or very
        % good) solution to the problem at hand. When set to logical true (1), offset is inserted in
        % the initial population. Thus, it is guaranteed that objective value at solution is not
        % greater than objective value at that, previously known, good point.
        if (options.trustoffset)
            X(1, :) = options.initoffset'; 
        end
    end
    
    % Initial velocities.
    % Velocities are initialized uniformly in [-VSPANINIT, VSPANINIT].
    if any(isnan(options.vspaninit))
        error('mrr:myoptim:pso:pso:optionserr:vspaninit', ...
                'VSPANINIT option must not contain NaN entries.');
    elseif isscalar(options.vspaninit)
        V = (rand(options.npart, nvars)-0.5)*2*options.vspaninit;
    else
        if (length(options.vspaninit) ~= size(options.vspaninit, 1)) || ...
           (length(options.vspaninit) ~= nvars)
            error('mrr:myoptim:pso:pso:optionserr:vspaninit', ...
                'VSPANINIT option must be either scalar or column-vector of length NVARS');
        end
        V = (rand(options.npart, nvars)-0.5)*2.*repmat(options.vspaninit', options.npart, 1);
    end
      
    % Initial scores (objective values).
    % Initialization of the best personal score and position, as well as global best score and
    % position.
    Y = calcobjfunc(objfunc, X);
    Ybest = Y;                      % The best individual score for each particle - initialization.
    Xbest = X;                      % The best individual position for each particle - 
                                    % initialization.
    [GYbest, gbest] = min(Ybest);   % GYbest is the best score within the entire swarm.
                                    % gbest is the index of particle that achived YGbest.
    gbest = gbest(1);               % In case when more than one particle achieved the best
                                    % score, we choose the one with the lowest index as the
                                    % best one.
                                    
    % These variables are used in testing mode only.
    tolbreak = ~isnan(options.globalmin);
    foundglobal = 0;
    if tolbreak && ~isscalar(options.globalmin)
        error('mrr:myoptim:pso:pso:optionserr:globalmin', ...
            'globalmin option, if specified, option must be a scalar value equal to the global minimum of the objective function');
    end
    
    % Initialization of the OUTPUT structure.
    % The output structure is filled and initialized differently depending on the OUTPUT_LEVEL
    % options, or equivalently depending on the output_level variable.
    if output_level >= 0
        % NONE log level
        output.itersno = options.niter;
        if output_level >= 1
            % LOW log level
            output.gbest_array = NaN*ones(options.niter+1, 1);
            output.gmean_array = NaN*ones(options.niter+1, 1);
            output.gworst_array = NaN*ones(options.niter+1, 1);
            output.gbest_array(1) = GYbest;
            output.gmean_array(1) = mean(Ybest);
            output.gworst_array(1) = max(Ybest);
            if output_level >= 2
                % MEDIUM log level
                output.gbestndx_array = NaN*ones(options.niter+1, 1);
                output.Xbest = NaN*ones(options.niter+1, nvars);
                output.gbestndx_array(1) = gbest;
                output.Xbest(1, :) = X(gbest, :);
                if output_level == 3
                    % HIGH log level
                    output.X = NaN*zeros(options.npart, nvars, options.niter+1);
                    output.X(:,:,1) = X;
                end
            end
        end
    end
    
    if options.verbose_period ~= 0
        disp 'PSO algorithm: Initiating the optimization process.'
    end

    % Denotes normal algorithm termination.
    exitflag = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The main loop.                                                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:options.niter
        
        % Verbosing, if neccessary.
        if options.verbose_period ~= 0
            if rem(iter, options.verbose_period) == 0
               disp(['iteration ', int2str(iter), '. best criteria = ', num2str(GYbest)]);
            end
        end
        
        % Calculating PSO parameters
        w = linrate(options.wf, options.wi, options.niter, 0, iter);
        cp = linrate(options.cbf, options.cbi, options.niter, 0, iter);
        cg = linrate(options.cgf, options.cgi, options.niter, 0, iter);

        % For later calculations only
        GXbest = repmat(Xbest(gbest, :), options.npart, 1);

        % Calculating speeds
        V = w*V + cp*rand(size(V)).*(Xbest-X) + cg*rand(size(V)).*(GXbest-X);
        V = min(vmax, abs(V)).*sign(V);

        % Population is moving
        X = X + V;
        Y = calcobjfunc(objfunc, X);

        % Calculating new individually best values
        mask = Y<Ybest;
        mask = repmat(mask, 1, nvars);
        Xbest = mask.*X +(~mask).*Xbest;
        Ybest = min(Y,Ybest);
        
        % Calculating new globally best value
        [GYbest, gbest] = min(Ybest);
        gbest = gbest(1);
        
        % Filling in the OUTPUT structure.
        if output_level >= 0
            % NONE log level
            if output_level >= 1
                % LOW log level
                output.gbest_array(iter+1)  = GYbest;
                output.gmean_array(iter+1)  = mean(Ybest);
                output.gworst_array(iter+1) = max(Ybest);
                if output_level >= 2
                    % MEDIUM log level
                    output.gbestndx_array(iter+1) = gbest;
                    output.Xbest(iter+1, :) = X(gbest, :);
                    if output_level == 3
                        % HIGH log level
                        output.X(:,:,iter+1) = X;
                    end
                end
            end
        end
        
        % The code used in testing mode only.
        if tolbreak && abs(GYbest - options.globalmin)<options.tol
            output.itersno = iter;
            foundglobal = 1;
            break
        end

    end
    
    if options.verbose_period ~= 0
        disp 'Optimization process finished.'
    end
    
    % Setting up the output variables.
    % X is set to be the final global best position, since that is the best position ever achieved
    % by any of the particles in the swarm. FVAL is set to be the value of the objective in X.
    x = Xbest(gbest, :); x = x(:);
    fval = GYbest;
    
    % The global moptimum has been found prior to achieving the maximal number of iteration.
    if foundglobal, exitflag = 1; 
    end
    
    % Plotting the algorithm behavior at each iteration.
    if options.plot
        r = 0:options.niter;
        figure
        plot(r, output.gbest_array, 'k.', r, output.gmean_array, 'r.', r, output.gworst_array, 'b.');
        str = sprintf('Best objective value : %g', fval);
        title(str);
        legend({'best objective', 'mean objective', 'worst objective'})
    end
    
end

function Y = calcobjfunc(func, X)
% CALCOBJFUNC   A helper function used to calculate objective function value for a series of points.
np = size(X,1);
Y = zeros(np,1);
for i = 1:np
    Y(i) = func(X(i,:));
end

function opts = getDefaultOptions
% GETDEFAULTOPTIONS     Returns a structure containing the default options.
%
%   This function, in fact, defines default values of the options within the options structure.
opts.npart          = 30;       % The number of particles.
opts.niter          = 100;      % The number of iterations.
opts.cbi            = 2.5;      % Initial value of the individual-best acceleration factor.
opts.cbf            = 0.5;      % Final value of the individual-best acceleration factor.
opts.cgi            = 0.5;      % Initial value of the global-best acceleration factor.
opts.cgf            = 2.5;      % Final value of the global-best acceleration factor.
opts.wi             = 0.9;      % Initial value of the inertia factor.
opts.wf             = 0.4;      % Final value of the inertia factor.
opts.vmax           = Inf;      % Absolute speed limit. It is the primary speed limit.
opts.vmaxscale      = NaN;      % Relative speed limit. Used only if absolute limit is unspecified.
opts.vspaninit      = 1;        % The initial velocity span. Initial velocities are initialized 
                                % uniformly in [-VSPANINIT, VSPANINIT]. 
opts.initoffset     = 0;        % Offset of the initial population.
opts.initspan       = 1;        % Span of the initial population.
opts.trustoffset    = 0;        % If set to 1 (true) and offset is vector, than the offset is 
                                % believed to be a good solution candidate, so it is included in 
                                % the initial swarm.
opts.initpopulation = NaN;      % The user-suplied initial population. If this is set to something
                                % meaningfull, then INITSPAN, INITOFFSET and TRUSTOFFSET are
                                % ignored.
opts.verbose_period = 10;       % The verbose period, i.e. the number of iterations after which the 
                                % results are prompted to the user. If set to 0, then verbosing is
                                % skipped.
opts.plot           = 0;        % If set to 1, evolution of the gbest is ploted to the user after
                                % the optimization process. The objective value of the best, mean
                                % and worse particle troughout the optimization process are plotted
                                % in the single graph.
opts.output_level   = 'low';    % The output log level. Possible values are: 'none', 'low', 
                                % 'medium', 'high'.
opts.globalmin      = NaN;      % Global minimum, used for testing only
opts.tol            = 1e-6;     % Precision tolerance, used for testing only

function x = linrate(xmax, xmin, tmax, tmin, t)
% LINRATE   Linear interpolation of value X in instant T, defined by previously known points
%           (tmin, xmin), (tmax, xmax)
x = xmin + ((xmax-xmin)/(tmax-tmin))*(tmax-t);
