function k = keipy_read_nc(filename) %,yr_temp)

    % ugly code assumes "*.nc" file is coming in here...
    matfilename = filename;
    matfilename(end-1:end+1) = 'mat';

    if (exist(matfilename,'file'))
        read_mat = 1;
%         if (exist(filename,'file'))
%             mf = dir(matfilename);
%             nf = dir(filename);
%             if (nf.datenum > mf.datenum)
%                 read_mat = 0;
%             end
%         end

        if (read_mat)
            load(matfilename)

            % take me out after I've updated moorings_8 MAT files
            %k.flx = ncread(filename,'flx');
            %save(matfilename,'k');
            %days = round(length(k.atm_flux_to_ice_surface)/24)-1;
            %for i=1:days
            %    k.total_ice_freeze_daily(i) = sum(k.total_ice_freeze(i*24-23:i*24));
            %end

            % print oct 1 ice thickness  oct1
%             k.dn = kei_dn(k,yr_temp);
%             dn_oct1 = datenum(yr_temp,10,1);
%             idx = find(k.dn == dn_oct1);
%             iceh = sprintf('Ice thickness Oct 1: %f',k.hi(idx(1)));
%             disp(iceh);

            return
        end
    end


    recs = length(ncread(filename,'f_time'));
    drange = length(ncread(filename,'zm'));

    sal_ref = 34.6;

    rv = get_read_vars(drange,recs);

    % preallocate
    %for i = 1:length(rv)
    %    command = sprintf('k.%s(%i,%i) = 0;',rv{i}{1},rv{i}{2},rv{i}{3});
        %disp(command);
    %    eval(command);
    %end

    %read data
    for i = 1:length(rv)
        [ncid] =netcdf.open(filename,'NOWRITE');
        command=sprintf('[varid]=netcdf.inqVarID(ncid,''%s\'');',rv{i}{1});
        %fprintf(1,'reading %s...\n',rv{i}{1});
        eval(command);
        netcdf.close(ncid);
        if (varid >= 0)
            command=sprintf('k.%s=ncread(filename,''%s\'');',rv{i}{1},rv{i}{1});
            %fprintf(1,'reading %s...\n',rv{i}{1});
        else
            command=sprintf('k.%s=NaN;',rv{i}{1},rv{i}{1});
        end
        eval(command);
        % transpose keipy netcdf vars
        command = sprintf('test_dims = size(k.%s);',rv{i}{1});
        eval(command);
        if (test_dims(2) > 1)
            command=sprintf('k.%s = k.%s'';',rv{i}{1},rv{i}{1});
            eval(command);
        end

        clear functions

    end

    % plot commands expect these vars
    k.zgrid = k.zm;
    k.hsn = k.hs;
    k.step = [0:recs-1];

    %pH
    depth_dbar = repmat(k.zm*(-1),1,recs);
    [ix,iy] = size(k.S);
    k.pH = zeros(ix,iy)*NaN;
    for i=1:ix
        for j=1:iy
            % pH_total is not array-safe it seems
            k.pH(i,j) = pH_total(k.S(i,j),k.T(i,j),depth_dbar(i,j),k.PO4(i,j),k.SiO3(i,j),k.ALK(i,j),k.DIC(i,j));
        end
    end

    % pull interesting fluxes

    %k.deep_ice_production = k.focn(8,:);  % this is included in other terms below
    k.surface_ice_production = -k.focn(7,:)*334000;  %(convert kg to W)
    k.total_ice_production = k.fio(7,:)*334000;  % < 0 (W)
    % behavior of ice on ocn in response to everything - usually negative
    % b/c melting.  Positive if melted completely and didn't use up ML heat
    % actual melt flux
    k.ice_to_ocn_flux = k.fio(4,:);  % (total, already scaled by ice fraction)
    % effectively the mixed layer melt potential
    k.ocn_to_ice_flux = k.fio(8,:); % not scaled by ice fraction

    % freeze potential
    k.ocn_freeze_potenial = k.fio(7,:)*334000; % kg -> Watts

    k.atm_flux_to_ice_surface=double(ncread(filename,'atm_flux_to_ice_surface'));  % already scaled by fraction ice
    k.atm_flux_to_ocn_surface=double(ncread(filename,'atm_flux_to_ocn_surface'));  % already scaled by fraction ocean

    k.total_ice_freeze=double(ncread(filename,'total_ice_freeze'));
    k.total_ice_melt=double(ncread(filename,'total_ice_melt'));
    k.frazil_ice_volume=double(ncread(filename,'frazil_ice_volume'));
    k.congelation_ice_volume=double(ncread(filename,'congelation_ice_volume'));
    k.snow_ice_volume=double(ncread(filename,'snow_ice_volume'));
    k.snow_precip_mass=double(ncread(filename,'snow_precip_mass'));
    k.ice_ocean_bottom_flux_potential=double(ncread(filename,'ice_ocean_bottom_flux_potential'));
    k.ice_ocean_bottom_flux=double(ncread(filename,'ice_ocean_bottom_flux'));

    pressure = sw_pres(-k.zm,-67);


    % calc extra bulk properties
    [~, n_steps] = size(k.T);
    for i=1:n_steps
        [~,k1] = min(k.T(:,i));
        k1=max([70,k1]);...
        k.deep_pyc(i) = mixedl(k.T(k1:end,i),k.T(k1,i),-1.*(k1-1:401-0.5));
        dp_i = round(abs(k.deep_pyc(i)));  % only works cause 1m grid!
        k.deep_pyc_end(i) = upside_down_mixedl(k.T(:,i),k.T(drange-1,i),-1.*(0:drange-0.5));
        k.deep_pyc_slope(i) = (k.T(dp_i+20,i) - k.T(dp_i,i))/20;
        k.deep_pyc_t(i) = mean(k.T(dp_i:dp_i+20,i));
        k.cp(:,i) = CPSW(k.S(:,i),k.T(:,i),-k.zm);
        k.sigma(:,i) = sw_dens(k.S(:,i),k.T(:,i),pressure);
        Tf = sw_fp(k.S(:,i),pressure);
        ml_idx = floor(-k.zml(i));
        k.heat(:,i) = k.cp(:,i).*k.sigma(:,i).*(k.T(:,i)-Tf);
        k.heat_60(i) = sum(k.heat(1:60,i));
        k.heat_ml(i) = sum(k.heat(1:ml_idx,i));
        k.heat_gr_60(i) = sum(k.heat(61:drange-1,i));
        sal_ref_i = find(k.S(60:end,i) > sal_ref,1) + 59;
        if isempty(sal_ref_i)
            sal_ref_i = drange-1;
        end
        kc = (sal_ref - k.S(1:sal_ref_i,i))./sal_ref;
        k.fresh(i) = sum(kc);
        kc = (sal_ref - k.S(1:ml_idx,i))./sal_ref;
        k.fresh_ml(i) = sum(kc);

        k.cp_prev(:,i) = CPSW(k.Sprev(:,i),k.Tprev(:,i),-k.zm);
        sigma(:,i) = sw_dens(k.Sprev(:,i),k.Tprev(:,i),pressure);
        Tf = sw_fp(k.Sprev(:,i),pressure);
        k.heat_prev(:,i) = k.cp_prev(:,i).*sigma(:,i).*(k.Tprev(:,i)-Tf);

        if (i > 1)
            %boundarylayer = round(k.hmx(i));
            %boundarylayer_1 = round(k.hmx(i-1));
            %k.pbl_flux(i) = sum(k.heat(1:boundarylayer,i) - k.heat_prev(1:boundarylayer,i))/(60*60); %dheat/timestep = flux
            k.pbl_flux(i) = (sum(k.heat(1:60,i)) - sum(k.heat(1:60,i-1)))/(60*60); %dheat/timestep = flux
            k.pyc_flux(i) = k.pbl_flux(i) - k.atm_flux_to_ocn_surface(i) - k.ice_to_ocn_flux(i) + k.total_ice_production(i);
        end
    end


    days = round(length(k.atm_flux_to_ice_surface)/24)-1;
    for i=1:days
        k.fice_daily(i) = mean(k.fice(i*24-23:i*24));
        k.fice_mean24(i*24-23:i*24) = k.fice_daily(i);
        k.surface_ice_production_daily(i) = mean(k.surface_ice_production(i*24-23:i*24));
        k.surface_ice_production_mean24(i*24-23:i*24) = k.surface_ice_production_daily(i);
        k.total_ice_production_daily(i) = mean(k.total_ice_production(i*24-23:i*24));
        k.total_ice_production_mean24(i*24-23:i*24) = k.total_ice_production_daily(i);
        k.ice_to_ocn_flux_daily(i) = mean(k.ice_to_ocn_flux(i*24-23:i*24));
        k.ice_to_ocn_flux_mean24(i*24-23:i*24) = k.ice_to_ocn_flux_daily(i);
        k.ocn_to_ice_flux_daily(i) = mean(k.ocn_to_ice_flux(i*24-23:i*24));
        k.ocn_to_ice_flux_mean24(i*24-23:i*24) = k.ocn_to_ice_flux_daily(i);
        k.atm_flux_to_ice_surface_daily(i) = mean(k.atm_flux_to_ice_surface(i*24-23:i*24));
        k.atm_flux_to_ice_surface_mean24(i*24-23:i*24) = k.atm_flux_to_ice_surface_daily(i);
        k.atm_flux_to_ocn_surface_daily(i) = mean(k.atm_flux_to_ocn_surface(i*24-23:i*24));
        k.atm_flux_to_ocn_surface_mean24(i*24-23:i*24) = k.atm_flux_to_ocn_surface_daily(i);

        k.ocn_freeze_potenial_daily(i) = mean(k.ocn_freeze_potenial(i*24-23:i*24));
        k.ocn_freeze_potenial_mean24(i*24-23:i*24) = k.ocn_freeze_potenial_daily(i);

        k.ice_ocean_bottom_flux_daily(i)=mean(k.ice_ocean_bottom_flux(i*24-23:i*24));

        k.hmx_daily(i) = mean(k.hmx(i*24-23:i*24));
        k.zml_daily(i) = mean(k.zml(i*24-23:i*24));
        k.hi_daily(i) = mean(k.hi(i*24-23:i*24));
        k.hs_daily(i) = mean(k.hs(i*24-23:i*24));

        k.deep_pyc_daily(i) = mean(k.deep_pyc(i*24-23:i*24));
        k.deep_pyc_end_daily(i) = mean(k.deep_pyc_end(i*24-23:i*24));
        k.deep_pyc_slope_daily(i) = mean(k.deep_pyc_slope(i*24-23:i*24));
        k.deep_pyc_t_daily(i) = mean(k.deep_pyc_t(i*24-23:i*24));

        k.fresh_daily(i) = mean(k.fresh(i*24-23:i*24));
        k.fresh_ml_daily(i) = mean(k.fresh_ml(i*24-23:i*24));

        k.heat_60_daily(i) = mean(k.heat_60(i*24-23:i*24));
        k.heat_ml_daily(i) = mean(k.heat_ml(i*24-23:i*24));
        k.heat_gr_60_daily(i) = mean(k.heat_gr_60(i*24-23:i*24));

        k.total_ice_freeze_daily(i) = sum(k.total_ice_freeze(i*24-23:i*24));

        k.frazil_ice_volume_daily(i)=sum(k.frazil_ice_volume(i*24-23:i*24));
        k.congelation_ice_volume_daily(i)=sum(k.congelation_ice_volume(i*24-23:i*24));
        k.snow_ice_volume_daily(i)=sum(k.snow_ice_volume(i*24-23:i*24));
        k.total_ice_volume_daily(i)=k.frazil_ice_volume_daily(i) + k.congelation_ice_volume_daily(i) + k.snow_ice_volume_daily(i);

        k.turb_daily_150(i)=sum(sum(k.wT(1:150,i*24-23:i*24)));


        k.pbl_flux_daily(i) = mean(k.pbl_flux(i*24-23:i*24));
        k.pyc_flux_daily(i) = mean(k.pyc_flux(i*24-23:i*24));

        % find zero-degree depth
        k.T_daily(:,i) = mean(k.T(:,i*24-23:i*24),2);
        k.d0(i) = -1;
        j = 200;
        while (k.d0(i) == -1 ) && (j > 0)
            if k.T_daily(j,i) < 0
                k.d0(i) = j;
            end
            j=j-1;
        end

        if j == 0
            %disp( 'zero degree water not found in kei_read_plotting!!!')
            %disp(sprintf('filename: %s, day: %i',filename,i))
            [~,k.d0(i)] = min(k.T_daily(50:end,i));
            k.d0(i) = k.d0(i)+49;
            %return
        end

    end

    k.d0 = k.d0*(-1);

    % save for later
    save(matfilename,'k')


    return

end


function rv = get_read_vars(drange,recs)

rv = { ...
{'zm',drange,1}, ...
{'day',recs,1}, ...
{'hour',recs,1}, ...
{'wT',drange,recs}, ...
{'wS',drange,recs}, ...
{'T',drange,recs}, ...
{'S',drange,recs}, ...
{'hmx',recs,1}, ...
{'zml',recs,1}, ...
{'km',drange,recs}, ...
{'ks',drange,recs}, ...
{'fatm',9,recs}, ...
{'fao',9,recs}, ...
{'fai',9,recs}, ...
{'fio',9,recs}, ...
{'focn',9,recs}, ...
{'hi',recs,1}, ...
{'hs',recs,1}, ...
{'fice',recs,1}, ...
{'ns',recs,1}, ...
{'ni',recs,1}, ...
{'dzi',42,recs}, ...
{'dzs',26,recs}, ...
{'Ts',26,recs}, ...
{'Ti',42,recs}, ...
{'Si',42,recs}, ...
{'ks',drange,recs}, ...
{'kt',drange,recs}, ...
{'tot_prod',drange,recs}, ...
{'diatC',drange,recs}, ...
{'diatChl',drange,recs}, ...
{'spC',drange,recs}, ...
{'spChl',drange,recs}, ...
{'ALK',drange,recs}, ...
{'DIC',drange,recs}, ...
{'DOC',drange,recs}, ...
{'DOFe',drange,recs}, ...
{'DON',drange,recs}, ...
{'DOP',drange,recs}, ...
{'Fe',drange,recs}, ...
{'NH4',drange,recs}, ...
{'NO3',drange,recs}, ...
{'O2',drange,recs}, ...
{'PO4',drange,recs}, ...
{'SiO3',drange,recs}, ...
{'diatFe',drange,recs}, ...
{'spCaCO3',drange,recs}, ...
{'spFe',drange,recs}, ...
{'diazC',drange,recs}, ...
{'diazChl',drange,recs}, ...
{'diazFe',drange,recs}, ...
{'diatSi',drange,recs}, ...
{'diat_Fe_lim',drange,recs}, ...
{'diat_light_lim',drange,recs}, ...
{'graze_diat',drange,recs}, ...
{'graze_tot',drange,recs}, ...
{'zooC',drange,recs}, ...
{'Tprev',drange,recs}, ...
{'Sprev',drange,recs}, ...
{'date',recs,1}, ... % now reading forcing for same netcdf
{'tau_x',recs,1}, ...
{'tau_y',recs,1}, ...
{'qswins',recs,1}, ...
{'qlwdwn',recs,1}, ...
{'tz',recs,1}, ...
{'qz',recs,1}, ...
{'prain',recs,1}, ...
{'psnow',recs,1}, ...
{'msl',recs,1}, ...
{'h',recs,1}, ...
{'dustf',recs,1}, ...
{'divu',recs,1}, ...
{'ic',recs,1}, ...
{'ain',recs,1}, ...
{'aout',recs,1}
};


end


function my_CP = CPSW(S,T1,P0)

% ******************************************************************
% UNITS:
%       PRESSURE        P0       DECIBARS
%       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
%       SALINITY        S        (IPSS-78)
%       SPECIFIC HEAT   CPSW     J/(KG DEG C)
% ***
% REF: MILLERO ET AL,1973,JGR,78,4499-4507
%       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
% PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
% DEVELOPED BY FOFONOFF 1980.
% ***
% CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
% T = 40 DEG C, P0= 10000 DECIBARS

		% inputs
    %real :: S,T1,P0

    % local
    %real :: T,P,SR,A,B,C,CP0,CP1,CP2

%   check that temperature is above -2
    T = T1;
    T(T < -2.0) = -2.0;

%   SCALE PRESSURE TO BARS
    P=P0./10.0;
% ***
% SQRT SALINITY FOR FRACTIONAL TERMS
    SR = sqrt(abs(S));
% SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
    A = (-1.38385E-3.*T+0.1072763).*T-7.643575;
    B = (5.148E-5.*T-4.07718E-3).*T+0.1770383;
    C = (((2.093236E-5.*T-2.654387E-3).*T+0.1412855).*T -3.720283).*T+4217.4;
    CP0 = (B.*SR + A).*S + C;
% CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8.*T+2.0357E-6).*T-3.13885E-4).*T+1.45747E-2).*T -0.49592;
    B = (((2.2956E-11.*T-4.0027E-9).*T+2.87533E-7).*T-1.08645E-5).*T +2.4931E-4;
    C = ((6.136E-13.*T-6.5637E-11).*T+2.6380E-9).*T-5.422E-8;
    CP1 = ((C.*P+B).*P+A).*P;
% CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10.*T+2.5941E-8).*T+9.802E-7).*T-1.28315E-4).*T +4.9247E-3;
    B = (3.122E-8.*T-1.517E-6).*T-1.2331E-4;
    A = (A+B.*SR).*S;
    B = ((1.8448E-11.*T-2.3905E-9).*T+1.17054E-7).*T-2.9558E-6;
    B = (B+9.971E-8.*SR).*S;
    C = (3.513E-13.*T-1.7682E-11).*T+5.540E-10;
    C = (C-1.4300E-12.*T.*SR).*S;
    CP2 = ((C.*P+B).*P+A).*P;
% SPECIFIC HEAT RETURN
    my_CP = CP0 + CP1 + CP2;
    return
end

function Tf = TfrzC(S,Db)

%     Freezing point of water in degrees C at salinity S in PSU
%                                         and pressure Db in decibars
    %TfrzC = (-0.0575 +1.710523e-3 *sqrt(S) -2.154996e-4 *S) *S &
    %- 7.53e-4 *Db

	Tf = -0.054.*S - 7.53e-4 .*Db;


    return
end



