% Input validation function for adcp orientation method for use with an input parser
%
% p = inputParser;
% addRequired(p,'adcp_orientation_method',@check_adcp_orientation_method)

function isvalid = check_adcp_orientation_method(x)
    isvalid = false;
    valid_methods = check_proc_methods([],'all','adcp_orientation');
    msg = sprintf('\nValid ADCP orientation methods are:');
    for i = 1:length(valid_methods)
        msg = sprintf('%s\n  - %s',msg,valid_methods{i});
    end

    if ~isstr(x) && ~iscell(x)
        error(['ADCP orientation method must be a string or a cell array containing a string and numeric parameter.' ...
               msg]);
    elseif iscell(x)
        method_name = x{1};
        if ~ismember(method_name,valid_methods)
            error(sprintf('Method "%s" is not recognized.%s',method_name,msg))
        elseif contains(lower(method_name),'offset')
            if length(x) < 2 || ~isnumeric(x{2})
                error(sprintf('ADCP orientation method "%s" requires numerical offset.',method_name))
            end
        end
        isvalid = true;
    elseif isstr(x)
        method_name = x;
        if ~ismember(method_name,valid_methods)
            error(sprintf('Method "%s" is not recognized.%s',method_name,msg))
        elseif contains(lower(method_name),'offset')
            error(sprintf('ADCP orientation method "%s" requires numerical offset.',method_name))
        end
        isvalid = true;
    end
end
