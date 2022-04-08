classdef Vessel
    % Base class for making sections

    properties
        instruments
    end

    methods
        function obj = Vessel(instrument)
            % Constructor
            arguments (Repeating)
                instrument
            end

            % cell array of instruments
            obj.instruments = instrument;
        end

        function output = get_all_data(obj,start_stop)

          arguments
	    obj
            start_stop datetime
          end

          % if start_stop is empty return an empty structure
          if isempty(start_stop); output = struct([]); return; end

	  start = start_stop(1);
	  stop = start_stop(2);

          % loop through the instruments
          for ins = 1:length(obj.instruments)
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
          end
        end
      end
end
