classdef Tchain < Instrument
    % Tchain from deployment files
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


    methods
        function obj = Tchain(name,data_directory,variables)
            % Constructor
            arguments
                name string = ""
                data_directory = string.empty
                variables string = string.empty
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
                obj.deployments(i).start_time = min(deployment.dn);
                obj.deployments(i).end_time = max(deployment.dn);
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
                   compare_vars = any(variables == string(who(deployment)));

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
               dn = deployment.dn;
               idx = dn >= datenum(start) & ...
                   dn <= datenum(stop);
               nidx = find(idx);
               clear dn

               % loop through variables saving each to the output structure
                for var = variables

                    switch var

                        case "info"
                            data{i}.info = deployment.info;

                        case "pos"
                            data{i}.pos = deployment.pos;
                        
                        otherwise
                            data{i}.(var) = deployment.(char(var))(:,nidx);
                            
                    end  
                    
                end
            end

        end

    end
    
end
