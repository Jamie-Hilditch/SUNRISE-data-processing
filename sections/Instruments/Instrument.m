classdef (Abstract) Instrument
    % Abstract base class for SUNRISE instruments
    %
    % Constructor Arguments:
    %   name: string (default="") - the name of the instrument
    %   variables: string (default=string.empty) - string array containing
    %                           the name of the instrument variables to be
    %                           loaded. If empty all variables should be
    %                           loaded.
    %
    % Properties:
    %   name: string - see constructor
    %   variables: string - see constructor
    %
    % Abstract Properties:
    %   data_source: e.g. file , list of files, or directory
    %
    %
    % Abstract Methods:
    %   get_data(obj,start,stop): Get data struct or cell of structs
    %                             containing the variables defined
    %                             by variables between two datetimes

    properties
        name string
        variables (1,:) string
    end

    properties (Abstract)
        data_source
    end

    methods
        function obj = Instrument(name,variables)
            % Constructor
            arguments
                name string = ""
                variables (1,:) string = string.empty
            end

            obj.name = name;
            obj.variables = variables;
        end
    end

    methods (Abstract)
        data = get_data(obj,start,stop)
    end
end
