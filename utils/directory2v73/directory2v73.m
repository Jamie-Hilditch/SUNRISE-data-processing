function directory2v73(directory)
  % Convert all mat files in a directory to v7.3

  files = dir(fullfile(directory,"*.mat"));
  for i = 1:length(files)
      clear adcp
      fname = fullfile(files(i).folder,files(i).name);
      data = load(fname);
      save(fname,'-struct','data','-v7.3');
      disp(['Saved ' fname])

  end
end
