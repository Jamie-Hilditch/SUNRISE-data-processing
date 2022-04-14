classdef Vessel
    % Base class for making sections

    properties
        instruments
    end

    methods
        function obj = Vessel(instrument)
            % Constructor
            arguments (Repeating)
                instrument {mustBeInstrument}
            end

            % cell array of instruments
            obj.instruments = instrument;
        end

        function obj = add_instrument(obj,instrument)

          arguments
            obj
          end

          arguments (Repeating)
            instrument {mustBeInstrument}
          end

          obj.instruments = [obj.instruments, instrument];
        end

        function output = get_all_data(obj,start,stop)

          arguments
	          obj
            start datetime
            stop datetime
          end

	  output = struct;

          % loop through the instruments
          for ins = 1:length(obj.instruments)
            try
              instrument = obj.instruments{ins};
              % get data
              data = instrument.get_data(start,stop);
              % if empty then there was no data for this instrument in this section
              if isempty(data); continue; end
              % some instruments return structs and some return cells arrays of structs
              % structures go straight into the output
              if ~iscell(data)
                output.(instrument.name) = data;
                continue
              end
              % for the cell arrays we need to extract the structures
              % we should only have one deployment per section but let's check
              if length(data) == 1
                output.(instrument.name) = data{1};
              else
                for i = 1:length(data)
                  name = append(instrument.name,sprintf('_%02d',i));
                  output.(name) = data{i};
                end
              end
            catch ME
              warning(sprintf('MyCode:CaughtError:%s',ME.identifier), ...
              sprintf('Caught error while getting data from %s\nError message:\n%s', ...
              instrument.name,ME.message))
            end

          end
        end
      end
end

function mustBeInstrument(input)
  if ~isa(input,'Instrument')
    error('MyCode:NotAnInstrument', ...
    'Input arguments must be implementations (subclasses) of Instrument')
  end
end
