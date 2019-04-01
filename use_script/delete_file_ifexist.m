function f = delete_file_ifexist(filename)
% This file delete an external file if it exist. Return 0 if operation is sucessful, or 1 if it fails

  if exist(filename, 'file')==2
    delete(filename);
    f= 0;
  else
    f = 1;
  end
end
