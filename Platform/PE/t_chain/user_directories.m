function d = user_directories(cruise)

switch cruise
  case 'SUNRISE'
      
%     d = fullfile(getenv('DATA'), 'SUNRISE', 'Tchain');
    d = fullfile('/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/', 'SUNRISE', 'Tchain');
end
