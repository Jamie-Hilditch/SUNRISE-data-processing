classdef VMP_Combo < Instrument
    % VMP from combo file
    % The combined file contains all the VMP casts from a single ship
    %
    % Constructor Arguments
    %   name: string (default="") - the name of the instrument
    %   data_file: (default=string.empty) - filepath to the vmp combo file
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
        function obj = VMP_Combo(name,data_file,variables)
            %Constructor
            arguments
                name string = ""
                data_file = string.empty
                variables string = string.empty
            end

            obj = obj@Instrument(name,variables);
            obj.data_source = data_file;
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

            %%%%%%%%%%
            % FIX ME %
            %%%%%%%%%%
            % Save combo file with -struct and then load variables
            % individually
            % load in data
            vmp_data = load(obj.data_source).vmp_combo;

            % obj.variables defines which variables are to be included
            % if obj.variables is empty (default) we include all the variables
            if isempty(obj.variables)
               variables = string(fieldnames(vmp_data))';
            else
               variables = obj.variables;
            end

            % Construct the time index for the section
            dn = vmp_data.dn;
            idx = dn >= datenum(start) & dn <= datenum(stop);

            % Return if no datapoints
            if sum(idx) == 0; data = struct([]); return; end

            % loop through variables saving each to the output structure
            for var = variables

                % Throw a warning if the variable does not exist but continue
                if ~isfield(vmp_data,var)
                    warning('Variable %s not found in %s',var,obj.name);
                    continue
                end

                % All vmp combo variables are two dimensional with the
                % dimension as the second dimension except depth and
                % u_star_cint
                switch var

                    case "depth"
                        data.depth = vmp_data.depth;

                    case "u_star_cint"
                        data.u_star_cint = vmp_data.u_star_cint(idx,:);

                    otherwise
                        data.(var) = vmp_data.(var)(:,idx);
                end
            end

            % print a message
            fprintf('    ├── %s\n',obj.name)     
        end
    end
end
