classdef ADCP_Rhib < Instrument
    % Rhib ADCPs from deployment files
    % The combined file contains all the VMP and  CTD casts from a single ship
    %
    % Constructor Arguments
    %   name: string (default="") - the name of the instrument
    %   data_directory: (default=string.empty) - filepath to the
    %                                            deployments directory
    %   variables: string (default=string.empty) - string array containing 
    %                           the name of the instrument variables to be
    %                           loaded. If empty all variables should be
    %                           loaded.
    %
    % Properties:
    %   name: string - see constructor
    %   variables: string - see constructor
    %   data_source: string - data_directory from constructor
    %   deployments: struct - filepath, start_time, and end_time for each
    %                         deployment file
    %
    %
    % Methods:
    %   get_data(obj,start,stop): Get data struct containing fields defined
    %                             in variables between two datetimes 

    properties
        data_source
    end

    properties (SetAccess = protected)
        deployments = struct([])
    end

    properties (Access = protected)
        % seperate the variables by different types of indexing
        % variable_aliases are defined by static methods
        variable_aliases = {"u","v","w","depth","z","dn"} %#ok<*CLARRSTR> 
        no_time_dims = {"config","cell_depth"}
        time_other = {"ahrs_gyro","computed_heading","vessel_u","vessel_v"}
        other_time = {"time","heading","pitch","roll","nuc_time","bt_range",...
                        "bt_range","bt_vel","bt_time","vessel_heading",...
                        "vessel_lat","vessel_lon"}
        other_other_time = {"vel","echo_intens","corr","ahrs_rotm"}
    end    

    methods
        function obj = ADCP_Rhib(name,data_directory,variables)
            % Constructor
            arguments
                name string {mustBeTextScalar} = ""
                data_directory {mustBeFolder} = string.empty
                variables string {mustBeText} = string.empty
            end

            obj = obj@Instrument(name,variables);
            obj.data_source = data_directory;

            % Return if data_directory has not been defined
            if isempty(obj.data_source); return; end

            % Now we save the start and end time for each deployment
            files = dir(fullfile(obj.data_source,"*.mat"));
            for i = length(files):-1:1
                obj.deployments(i).path = fullfile(files(i).folder,files(i).name);
                deployment = matfile(obj.deployments(i).path);
                obj.deployments(i).start_time = min(deployment.time);
                obj.deployments(i).end_time = max(deployment.time);
                clear deployment
            end

        end

        function data = get_data(obj,start,stop)
            % get data between start and stop
            % Retrieve data from the data_file between the start and stop 
            % times. The variables included are defined by the variables
            % property and if this is empty all variables are included.
            %
            % Arguments
            %   start: datetime - start time of the section
            %   stop: datetime - end time of the section
            %
            % Output
            %   data: cell - cell array of structs, one per deployment,
            %           with fields defined by the variables property

            arguments
                obj
                start datetime
                stop datetime
            end
            
            % first find the correct deployment(s)
            deploy_idx = find([obj.deployments(:).start_time] <= datenum(stop) & ...
                [obj.deployments(:).end_time] >= datenum(start));


           % for each deployment return a struct, store each in a cell
           % array
           num_deploys = length(deploy_idx);
           data = cell(1,num_deploys);
           for i = 1:num_deploys

               % get matfile for this deployment
               deployment = matfile(obj.deployments(deploy_idx(i)).path);

               % get variables
               if isempty(obj.variables)
                   variables = string(who(deployment))';
               else
                   variables = obj.variables;

                   % check all these variables exist
                   all_variables = [string(who(deployment));obj.variable_aliases'];
                   compare_vars = any(variables == all_variables);

                   % throw a warning for variables that don't exist
                   unrecognized_vars = variables(~compare_vars);
                   for var = unrecognized_vars
                       warning('Variable %s not found in %s - file %s', ...
                           var,obj.name,deployment.Properties.Source);
                   end

                   % continue with variables that do exist
                   variables = variables(compare_vars);
               end

               % get time index
               time = deployment.time;
               idx = time >= datenum(start) & ...
                   time <= datenum(stop);
               nidx = find(idx);
               clear time

               % loop through variables saving each to the output structure
                for var = variables

                    switch var

                        % variable aliases are defined as static methods
                        % below
                        case obj.variable_aliases
                            data{i}.(var) = obj.(var)(deployment,nidx);

                        case obj.no_time_dims
                            data{i}.(var) = deployment.(var);

                        case obj.time_other
                            data{i}.(var) = deployment.(char(var))(nidx,:);

                        case obj.other_time
                            data{i}.(var) = deployment.(char(var))(:,nidx);

                        case obj.other_other_time
                            data{i}.(var) = deployment.(char(var))(:,:,nidx);
                        
                        otherwise
                            % unexpected variable
                            % may need to add to the static properties
                            % above
                            tmp = deployment.(var);
                            
                            % if not numeric save everything
                            if ~isnumeric(tmp); data{i}.(var) = tmp; continue; end

                            % if numeric try and find the time dimension
                            t_idx = find(size(tmp) == length(idx));
                            
                            % deal with the case of no time dimension
                            if isempty(t_idx);  data{i}.(var) = tmp; continue; end
                            
                            % index into time dimension
                            data{i}.(var) = subsasgn(tmp, ...
                                struct('type','()','subs',{ ...
                                [repmat({':'},1,t_idx(1)-1), ...
                                idx,...
                                repmat({':'},1,ndims(shifted)-t_idx(1))] ...
                                }),[]);
                            
                    end
                    
                end
            end

        end

    end
    
    methods (Static = true, Access = protected)
        function out = u(deployment,nidx)
            out = squeeze(deployment.vel(:,1,nidx));
        end

        function out = v(deployment,nidx)
            out = squeeze(deployment.vel(:,2,nidx));
        end

        function out = w(deployment,nidx)
            out = squeeze(deployment.vel(:,3,nidx));
        end

        function out = depth(deployment,~)
            out = deployment.cell_depth;
        end

        function out = z(deployment,~)
            out = -deployment.cell_depth;
        end

        function out = dn(deployment,nidx)
            out = deployment.nuc_time(:,nidx);
        end
    end
end