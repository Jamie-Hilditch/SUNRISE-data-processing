classdef ADCP_Ship_Combo < Instrument
    % Shipboard ADCP from combined file
    % The combined file contains all the ADCPs from a single ship and
    % therefore we require the ADCP name in addition to the file
    %
    % Constructor Arguments
    %   name: string (default="") - the name of the instrument
    %   data_file: (default=string.empty) - filepath
    %   variables: string (default=string.empty) - string array containing
    %                           the name of the instrument variables to be
    %                           loaded. If empty all variables should be
    %                           loaded.
    %
    % Properties:
    %   name: string - see constructor
    %   variables: string - see constructor
    %   data_source: string - data_file from constructor
    %
    %
    % Methods:
    %   get_data(obj,start,stop): Get data struct containing fields defined
    %                             in variables between two datetimes

    properties
        data_source
    end

    methods
        function obj = ADCP_Ship_Combo(name,data_file,variables)
            % Constructor
            arguments
                name string = ""
                data_file = string.empty
                adcp_name string = ""
                variables string = string.empty
            end

            % check file and adcp_name are valid
            varlist = who(matfile(data_file));
            if ~isempty(varlist) && ~any(strcmp(varlist,adcp_name))
                cell_adcps = join(varlist,', ');
                error('MATLAB:mycode:variableNotFound', ...
                    ['Could not find %1s in file %2s.\n', ...
                    'Available ADCPs are %s\n'], ...
                    adcp_name,data_file,cell_adcps{:})
            end

            obj = obj@Instrument(name,variables);

            obj.data_source = data_file;
            obj.adcp_name = adcp_name;
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
            %   data: struct - struct with fields defined by the
            %           variables property

            arguments
                obj
                start datetime
                stop datetime
            end

            % load in data
            adcp_data = matfile(obj.data_source);;

            % obj.variables defines which variables are to be included
            % if obj.variables is empty (default) we include all the variables
            % get variables
            if isempty(obj.variables)
                variables = string(who(adcp_data))';
            else
                variables = obj.variables;

                % check all these variables exist
                compare_vars = any(variables == string(who(adcp_data)));

                % throw a warning for variables that don't exist
                unrecognized_vars = variables(~compare_vars);
                for var = unrecognized_vars
                    warning('Variable %s not found in %s - file %s', ...
                        var,obj.name,adcp_data.Properties.Source);
                end

                % continue with variables that do exist
                variables = variables(compare_vars);
            end

            % Construct the time index for the section
            dn = adcp_data.dn;
            idx = dn >= datenum(start) & dn <= datenum(stop);

            % Return if no datapoints
            if sum(idx) == 0; data = struct([]); return; end

            % loop through variables saving each to the output structure
            for var = variables

                if var == "dn"
                  data.dn = adcp_data.dn(idx,:);
                end

                % All other ship ADCP variables are two dimensional with the
                % dimension as the second dimension
                data.(var) = adcp_data.(var)(:,idx);

            end
            % print a message
            fprintf('    ├── %s\n',obj.name)
        end

    end
end
