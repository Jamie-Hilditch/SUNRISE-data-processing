classdef FTMET_mat < Instrument

  properties
    data_source
  end

  methods
    function obj = FTMET_mat(name,datafile,variables)
      %constructor
      arguments
        name = ""
        datafile = string.empty
        variables = string.empty
      end

      obj = obj@Instrument(name,variables);
      obj.data_source = datafile;
    end

    function data = get_data(obj,start,stop)

      arguments
        obj
        start datetime
        stop datetime
      end

      % open file
      ftmet = matfile(obj.data_source);

      % obj.variables defines which variables are to be included
      % if obj.variables is empty (default) we include all the variables
      % get variables
      if isempty(obj.variables)
          variables = string(who(ftmet))';
      else
          variables = obj.variables;

          % check all these variables exist
          compare_vars = any(variables == string(who(ftmet)));

          % throw a warning for variables that don't exist
          unrecognized_vars = variables(~compare_vars);
          for var = unrecognized_vars
              warning('Variable %s not found in %s - file %s', ...
                  var,obj.name,ftmet.Properties.Source);
          end

          % continue with variables that do exist
          variables = variables(compare_vars);
      end

      % Construct the time index for the section
      dn = ftmet.time;
      nidx = find(dn >= datenum(start) & dn <= datenum(stop));

      % loop through variables saving each to the output structure
      for var = variables

        % only one variable is not a vector
        if var == "readme"
          data.readme = ftmet.readme;
          continue
        end

        % index and save data in all the other cases
        data.(var) = ftmet.(char(var))(nidx,:);

      end

      % print a message
      fprintf('    ├── %s\n',obj.name)
    end
  end
end
