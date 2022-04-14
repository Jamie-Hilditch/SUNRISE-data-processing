classdef ADCP_Ship_Combo < Instrument
    % Shipboard ADCP from combined file
    % The combined file contains all the ADCPs from a single ship and
    % therefore we require the ADCP name in addition to the file
    %
    % Constructor Arguments
    %   name: string (default="") - the name of the instrument
    %   data_file: (default=string.empty) - filepath
    %   adcp_name: string (default="") - name of the adcp as defined in the
    %                                    combo file
    %   variables: string (default=string.empty) - string array containing
    %                           the name of the instrument variables to be
    %                           loaded. If empty all variables should be
    %                           loaded.
    %
    % Properties:
    %   name: string - see constructor
    %   variables: string - see constructor
    %   data_source: string - data_file from constructor
    %   adcp_name: string - see constructor
    %
    %
    % Methods:
    %   get_data(obj,start,stop): Get data struct containing fields defined
    %                             in variables between two datetimes

    properties
        data_source
        adcp_name string
    end

    methods
        function obj = ADCP_Ship_Combo(name,data_file,adcp_name,variables)
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
            adcp_data = load(obj.data_source,obj.adcp_name).(obj.adcp_name);

            % obj.variables defines which variables are to be included
            % if obj.variables is empty (default) we include all the variables
            if isempty(obj.variables)
               variables = string(fieldnames(adcp_data))';
            else
               variables = obj.variables;
            end

            % FIXME: All processed files should have a datenum field.
            % Workaround for now.
            if ~isfield(adcp_data,'dn')
                adcp_data.dn = datenum(adcp_data.time');
            end

            % Construct the time index for the section
            dn = adcp_data.dn;
            idx = dn >= datenum(start) & dn <= datenum(stop);

            % Return if no datapoints
            if sum(idx) == 0; data = struct([]); return; end

            % loop through variables saving each to the output structure
            for var = variables

                % Throw a warning if the variable does not exist but continue
                if ~isfield(adcp_data,var)
                    warning('Variable %s not found in %s',var, ...
                        append(obj.name, ' - ', obj.adcp_name));
                    continue
                end

                % All ship ADCP variables are two dimensional with the
                % dimension as the second dimension
                data.(var) = adcp_data.(var)(:,idx);

            end
            % print a message
            fprintf('    ├── %s\n',obj.name)
        end

    end
end
