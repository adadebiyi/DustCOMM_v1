function [f] = write_netcdf (filename, ndimension, dim_name, varname, vardata)
% =========
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% =========
% This function takes the filename, variables and its data, and write it in a netcdf. The idea is to make things easier.
% Return 1 if all goes well, and 0 is there is a problem.
% ndimension = is the number of dimension
% dim_name = is the name of the dimensions
% varname  = is the name of the variables to be stored
% vardata = data for the variable arrange in the order of the dimension name

% =========  

  if (exist('filename','var') && exist('ndimension','var') && exist('dim_name','var') && exist('varname','var') && exist('vardata','var') )
    ndim_size  = size(vardata);
    switch ndimension
    case 1 % dims array
        ndim_size  = length(vardata);
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name),ndim_size});
      case 2 % dim array
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name(1)),ndim_size(1),char(dim_name(2)),ndim_size(2)});
      case 3 % dim array
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name(1)),ndim_size(1),char(dim_name(2)),ndim_size(2),char(dim_name(3)),ndim_size(3)});
      case 4 % dim array
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name(1)),ndim_size(1),char(dim_name(2)),ndim_size(2),char(dim_name(3)),ndim_size(3),char(dim_name(4)),ndim_size(4)});
      case 5 % dim array
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name(1)),ndim_size(1),char(dim_name(2)),ndim_size(2),char(dim_name(3)),ndim_size(3),char(dim_name(4)),ndim_size(4),char(dim_name(5)),ndim_size(5)});
      case 6 % dim array
        nccreate(char(filename),char(varname),'Dimensions',{char(dim_name(1)),ndim_size(1),char(dim_name(2)),ndim_size(2),char(dim_name(3)),ndim_size(3),char(dim_name(4)),ndim_size(4),char(dim_name(5)),ndim_size(5),char(dim_name(6)),ndim_size(6)});
      end %switch, i

    ncwrite(char(filename),char(varname),vardata)
    f=1;
  else
    disp('Write_netcdf: not all dimensions are represented: size(dim_name,2)~=nsize')
    f=0;
  end
