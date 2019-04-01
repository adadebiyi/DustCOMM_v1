 function [model_column_loading,model_D_lower,model_D_upper,N] = scale_model_no_3d(load_column_mass,D_lower,D_upper,model_name,ref_fraction,Dlw_limit,Dup_limit)
   
% Author: Adeyemi Adebiyi
% Affliation: University of California Los Angeles

% Correcting each model to be between 0.2 and 20 microns.
% -----

% =====================
%  define some global variables
  define_global_variables (true);
  global gPSD_OBS_median D_OBS

% set the dimensions
  ndim = size(load_column_mass);
  nlon = ndim(2);
  nlat = ndim(3);
	nlev = ndim(4);

% =====================
% Now begin scaling each model
% =====================
   if any(strcmp(model_name, 'GISS'))
% % GISS
     disp(['Scaling model data for ... GISS']);

   % copy the  model diameter limits
     N = max(size(find(D_lower(:) ~= 0)));
     Dlower = D_lower(1:N);
     Dupper = D_upper(1:N);

   % ; scale last
     scale_factor = PSD_V_integral(gPSD_OBS_median,D_OBS,16,Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,8,16);
     load_1620 = load_column_mass(N-1,:,:,:)*scale_factor;

   % copy & store
     dust_column_burden_kg_m2 = squeeze(load_column_mass(1:N-1,:,:,:));
     dust_column_burden_kg_m2(N,:,:,:) = load_1620;

     model_D_lower(1:N) = [Dlower(1:N-1), Dupper(N-1)];
     model_D_upper(1:N) = [Dupper(1:N-1), Dup_limit];

   %Normalize and store
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;
     model_column_loading(5,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(5,:,:,:)./CM;
     model_column_loading(6,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(6,:,:,:)./CM;
     model_column_loading(7,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(7,:,:,:)./CM;
     model_column_loading(8,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(8,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_1620;
   end

   if any(strcmp(model_name, 'WRFChem'))
     % =====================
 % %  WRFChems
     disp('Scaling model data for ... WRFChems');

     %   copy the  model diameter limits
     N = max(size(find(D_lower(:) ~= 0)));
     Dlower = D_lower(1:N);
     Dupper = D_upper(1:N);

     % Scale last
     numerat = PSD_V_integral(gPSD_OBS_median,D_OBS,Dupper(end),Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(end),Dupper(end));
     denominat = PSD_V_integral(gPSD_OBS_median,D_OBS,8,16)/PSD_V_integral(gPSD_OBS_median,D_OBS,4,8);
     scale_factor = numerat/denominat;
     scale_last = scale_factor*(ref_fraction(end-1,:,:,:)./ref_fraction(end-2,:,:,:));
     load_last = squeeze(load_column_mass(N,:,:,:)).*squeeze(scale_last);

   % scale first
     scale_first = PSD_V_integral(gPSD_OBS_median,D_OBS,Dlw_limit,Dupper(1))/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(1),Dupper(1));
     load_first = load_column_mass(1,:,:,:)*scale_first;

   % copy & store
     N = N+1;
     model_D_lower(1:N) = [Dlw_limit,Dlower(2:end),Dupper(end)];
     model_D_upper(1:N) = [Dupper,Dup_limit];

     dust_column_burden_kg_m2(1,:,:,:) = squeeze(load_first);
     dust_column_burden_kg_m2(2:N-1,:,:,:) = load_column_mass(2:N-1,:,:,:);
     dust_column_burden_kg_m2(N,:,:,:) = load_last;

   % normalize
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;
     model_column_loading(5,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(5,:,:,:)./CM;
     model_column_loading(6,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(6,:,:,:)./CM;
     model_column_loading(7,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(7,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_first load_last;

   end

   if any(strcmp(model_name, 'CESM'))

   % =====================
   % % CESM
     disp('Scaling model data for ... CESM');

     %   copy the  model diameter limits
     N = max(size(find(D_lower(:) ~= 0)));
     Dlower = D_lower(1:N);
     Dupper = D_upper(1:N);

     %scale last
     numerat = PSD_V_integral(gPSD_OBS_median,D_OBS,Dupper(end),Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(end),Dupper(end));
     denominat = PSD_V_integral(gPSD_OBS_median,D_OBS,8,16)/PSD_V_integral(gPSD_OBS_median,D_OBS,4,8);
     scale_factor = numerat/denominat;
     scale_last = scale_factor*(ref_fraction(end-1,:,:,:)./ref_fraction(end-2,:,:,:));
     load_last = squeeze(load_column_mass(N,:,:,:)).*squeeze(scale_last);

   %Scale first
     scale_first = PSD_V_integral(gPSD_OBS_median,D_OBS,Dlw_limit,Dupper(1))/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(1),Dupper(1));
     load_first = load_column_mass(1,:,:,:)*scale_first;

     % copy and store
     N = N+1;
     model_D_lower(1:N) = [Dlw_limit,Dlower(2:end),Dupper(end)];
     model_D_upper(1:N) = [Dupper,Dup_limit];

     dust_column_burden_kg_m2(1,:,:,:) = squeeze(load_first);
     dust_column_burden_kg_m2(2:N-1,:,:,:) = load_column_mass(2:N-1,:,:,:);
     dust_column_burden_kg_m2(N,:,:,:) = load_last;

   %normalize
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;
     model_column_loading(5,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(5,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_first load_last;

   end

   if any(strcmp(model_name, 'GEOSChem'))
     % =====================
   % % GEOSChem
     disp('Scaling model data for ... GEOSChem');

     %   copy the  model diameter limits
     N = max(size(find(D_lower(:) ~= 0)));
     Dlower = D_lower(1:N);
     Dupper = D_upper(1:N);

   %scale last
     numerat = PSD_V_integral(gPSD_OBS_median,D_OBS,Dupper(end),Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(end),Dupper(end));
     denominat = PSD_V_integral(gPSD_OBS_median,D_OBS,8,16)/PSD_V_integral(gPSD_OBS_median,D_OBS,4,8);
     scale_factor = numerat/denominat;
     scale_last = scale_factor*(ref_fraction(end-1,:,:,:)./ref_fraction(end-2,:,:,:));
     load_last = squeeze(load_column_mass(N,:,:,:)).*squeeze(scale_last);

     %copy and store
     N = N + 1;
     dust_column_burden_kg_m2 = squeeze(load_column_mass(1:N-1,:,:,:));
     dust_column_burden_kg_m2(N,:,:,:) = load_last;

     model_D_lower(1:N) = [Dlower(1:N-1), Dupper(N-1)];
     model_D_upper(1:N) = [Dupper(1:N-1), Dup_limit];

     %normalize
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;
     model_column_loading(5,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(5,:,:,:)./CM;
     model_column_loading(6,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(6,:,:,:)./CM;
     model_column_loading(7,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(7,:,:,:)./CM;
     model_column_loading(8,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(8,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_last;

   end

   if any(strcmp(model_name, 'CNRM'))
     % =====================
   % % GEOSChem
     disp('Scaling model data for ... CNRM-ESM');

     %   copy the  model diameter limits
     model_N = max(size(D_lower(:) ~= 0));
     Dlower = D_lower(2:model_N-1);
     Dupper = D_upper(2:model_N-1);

% ; scale last
      scale_factor = PSD_V_integral(gPSD_OBS_median,D_OBS,Dupper(end),Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(end),Dupper(end));
      load_last = squeeze(load_column_mass(model_N-1,:,:,:)).*squeeze(scale_factor);

%copy and store
     N = max(size(Dupper))+1;
     dust_column_burden_kg_m2(1:N-1,:,:,:) = squeeze(load_column_mass(2:model_N-1,:,:,:));
     dust_column_burden_kg_m2(N,:,:,:) = load_last;

     model_D_lower(1:N) = [Dlower(1:N-1), Dupper(N-1)];
     model_D_upper(1:N) = [Dupper(1:N-1), Dup_limit];

     %normalize
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;
     model_column_loading(5,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(5,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_last;

   end


   if any(strcmp(model_name, 'IMPACT'))
     % =====================
     disp('Scaling model data for ... IMPACT');

%   copy the  model diameter limits
     Nstart = 2;
     model_N = max(size(D_lower(:) ~= 0));
     Dlower = D_lower(Nstart:model_N);
     Dupper = D_upper(Nstart:model_N);

% Scale last
     numerat = PSD_V_integral(gPSD_OBS_median,D_OBS,Dupper(end),Dup_limit)/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(end),Dupper(end));
     denominat = PSD_V_integral(gPSD_OBS_median,D_OBS,8,16)/PSD_V_integral(gPSD_OBS_median,D_OBS,4,8);
     scale_factor = numerat/denominat;
     scale_last = scale_factor*(ref_fraction(end-1,:,:,:)./ref_fraction(end-2,:,:,:));
     load_last = squeeze(load_column_mass(end,:,:,:)).*squeeze(scale_last);

    %Scale first
      scale_first = PSD_V_integral(gPSD_OBS_median,D_OBS,Dlw_limit,Dupper(1))/PSD_V_integral(gPSD_OBS_median,D_OBS,Dlower(1),Dupper(1));
      load_first = load_column_mass(Nstart,:,:,:)*scale_first;

     %copy and store
     N = max(size(Dupper))+1;
     dust_column_burden_kg_m2(1,:,:,:) = load_first;
     dust_column_burden_kg_m2(Nstart:N-1,:,:,:) = squeeze(load_column_mass(Nstart+1:end,:,:,:));
     dust_column_burden_kg_m2(N,:,:,:) = load_last;

     model_D_lower(1:N) = [Dlw_limit, Dlower(Nstart:end),Dupper(end)];
     model_D_upper(1:N) = [Dupper, Dup_limit];

     %normalize
     CM = sum(dust_column_burden_kg_m2,1,'omitnan');
     model_column_loading(1,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(1,:,:,:)./CM;
     model_column_loading(2,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(2,:,:,:)./CM;
     model_column_loading(3,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(3,:,:,:)./CM;
     model_column_loading(4,1:nlon,1:nlat,1:nlev) = dust_column_burden_kg_m2(4,:,:,:)./CM;

     clear dust_column_burden_kg_m2 CM load_last;
   end

 end
