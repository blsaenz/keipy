% plt_nc -- modifed into function, Ben Saenz 7/2011

% !!! hardwired for hourly time step !!!!

% quick ice-ocean plots

function [t2_crop,surf_chla,total_diatC,total_spC,total_zooC] = keipy_plot_ecoiceline_dn(kei,d,range,year)
    
	% hardwired for hourly time step
	time2=zeros(length(kei.step),1);
	dtday=1/24.; 
	stday=mod(kei.day(1),365);
	time2 = double(stday) + double(kei.step).*dtday;
    time2_dn=kei_dn(kei,year);
	
	xlen=length(kei.zml);
	x=[0:xlen-1];
    
    cb_expand = 0.07;
	
	fig = figure;
    if(~usejava('Desktop'))
       set(fig,'visible','off');
    end
    set(fig,'Position',[250 150 700 920]);
	%pause;
	
    datastart = range(1);
    dataend = range(2);
    drange = datastart:dataend;
    t2_crop = time2(drange);
    t2dn_crop = time2_dn(drange);
    [y_s,m_s,d_s,h_s,mi_s,s_s] = datevec(t2dn_crop(1));
    [y_e,m_e,d_e,h_e,mi_e,s_e] = datevec(t2dn_crop(end));
    tick_locations = datenum(y_s,[m_s+1:(y_e-y_s)*12+m_e],1);
    

    load seawifs_LTER_200_300
    %plot(dn,chla,'-om','MarkerSize',3); hold on;
    seawifs_range = find((dn >= t2dn_crop(1)) & (dn <= t2dn_crop(end)));
    seawifs_time2_crop_dn = dn(seawifs_range);
    seawifs_time2_crop = seawifs_time2_crop_dn - datenum(year,0,0);
    

    % prepare valid ice/snow data
    % =====================================================================
    ns=kei.ns;
    ni=kei.ni;
    icet=kei.Ti-273.15;
    snowt=kei.Ts-273.15;
    icez=kei.dzi;
    snowz=kei.dzs;
    imask=zeros(size(icet));
    smask=zeros(size(snowt));
    for i=1:xlen
       imask(1:ni(i),i) = 1;
       smask(1:ns(i),i) = 1;
    end
    icet=icet.*imask;
    icez=icez.*imask;
    ices=kei.Si.*imask;
    snowt=snowt.*smask;
    snowz=snowz.*smask;
    total_hi = kei.hi+kei.hsn;

    % =====================================================================
      
    thickness = icez;
    z_ice = icez;
    xmax=xlen;    
    
    [y_layers,xmax] = size(icez);
    y_total = y_layers*5;
    %y_plot = fix(0.9*y_total);
    y_plot = y_total;
    
    %max_depth = max(total_hi);
    max_depth = 1.2;
    y_scale = y_plot/max_depth;

    th_snow = snowz;
    max_snow_depth = max(kei.hsn);
    y_scale_snow = 0.9*y_plot/max_snow_depth;

    
    % =====================================================================    
    
	clf
    
    subplot(7,1,2)
	plot(t2dn_crop,-kei.hmx(drange),'k'); hold on;
	plot(t2dn_crop,kei.zml(drange),'k--')
    title('Mixed Layer Depth (m)')
    ylim([-d 0]);
    
    set(gca,'XTick',tick_locations)
    datetick('x','mmm', 'keepticks')

    subplot(7,1,3)
	plot(t2dn_crop,-kei.hi(drange),'k'); hold on;
	plot(t2dn_crop,kei.hsn(drange),'k'); hold on;
	plot(t2dn_crop,zeros([1, length(drange)]),'k--');
    title('Sea Ice and Snow Thickness (m)')
    ylim([-1.2 0.6]);
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
   
    
    %plot(t2_crop,kei.fice(drange),'c'); hold on;
    %plot(t2_crop(1:length(frange)),f.forcing.ic(frange),'c--');
    %title('Ice concentration');legend('Model Ice','SSM/I ice');
  
    
    subplot(7,1,5)
    wind = sqrt(kei.tau_x.^2+kei.tau_y.^2);
    plot(t2dn_crop,wind(drange),'m');hold on;
    title('Wind Speed (m/s)')
    %for i=1:fix(length(time2)/24)
    %    ii = i*24;
    %    time24(i) = time2(ii-11);
    %    wind24(i) = mean(wind(ii-23:ii));
    %end
    wind1=wind;
    for ii=1:length(t2dn_crop)/24-1
        time24(ii) = t2dn_crop(ii*24-12);
        wind24(ii) = mean(wind1(ii*24-23:ii*24));
    end
    plot(time24,wind24,'b');
    axx = axis;axx(4)=20;axis(axx);
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
    legend('Hourly','24 Mean');
    
	subplot(7,1,7);
    surf_chla = mean(kei.diatChl(1:2,:)+kei.spChl(1:2,:));
    total_diatC = sum(kei.diatC(1:d,:));
    total_spC = sum(kei.spC(1:d,:));
    total_zooC = sum(kei.zooC(1:d,:));
    
    plot(t2dn_crop,surf_chla(drange),'g');
    hold on;
    plot(seawifs_time2_crop_dn,chla(seawifs_range),'-om','MarkerSize',3);
    ylim([0 5]);
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
    title('Surface Chl a (mg m^-^2)'); legend('Model','SeaWIFS');
 	xlabel('Time (day)')
    
    %axx = findall(fig,'type','axes');
    axx = findobj(fig,'type','axes','-not','Tag','legend','-not','Tag','colorbar');
    for i=1:length(axx)
        xlim(axx(i),[min(t2dn_crop);max(t2dn_crop)]);
        %set(axx(i),'XTick',[365-forcestart/24;720-forcestart/24;1095-forcestart/24]);
        %set(axx(i),'XTicklabel',[365 730 1095]);
    end
    
    s1 = subplot(7,1,1);
    
    imagesc(t2dn_crop,-kei.zgrid(1:d),kei.T(1:d,drange),[-2 5]); hold on;
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
	ylabel('Depth')
	title('Temperature')

    s1pos = get(s1,'position');
	colorbar('EastOutside');
    s11pos = get(s1,'position');     
    set(s1,'position',s1pos);

    
    s1=subplot(7,1,6);
    imagesc(t2dn_crop,-kei.zgrid(1:d),kei.Fe(1:d,drange),[0 6e-4]);
	ylabel('Depth')
    title ('Fe (0-0.6 nM)')
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
    
    s1pos = get(s1,'position');
	colorbar('EastOutside');
    s11pos = get(s1,'position');     
    %s1pos(2) = s11pos(2) - cb_expand/2;
    %s1pos(4) = s11pos(4) + cb_expand;
    set(s1,'position',s1pos);

    s2 = subplot(7,1,4);
    y_ref = 1;
    image_data = zeros(y_total,xmax);
    for i=drange
        y_ref=1;
        for j=1:ni(i)
            y_ref_new = y_ref+fix(icez(j,i)*y_scale);
            image_data(y_ref:y_ref_new,i) = ices(j,i);
            y_ref=y_ref_new+1;
        end
    end
    [idy, idx] = size(image_data);
    
    ice_depth_y = (1:idy)./idy.*max_depth;
    
  	imagesc(t2dn_crop,ice_depth_y,image_data,[0 20]); %colorbar;
    %pcolor(flipud(image_data));caxis([0 20]);shading flat;colorbar;
    title('Ice Salinity (0 - 20 ppt)');
    set(gca,'XTick',tick_locations)
    datetick('x','mmm','keeplimits', 'keepticks')
    
    s2pos = get(s2,'position');
	colorbar('EastOutside');
    s11pos = get(s2,'position');     
    %s1pos(2) = s11pos(2) - cb_expand/2;
    %s1pos(4) = s11pos(4) + cb_expand;
    set(s2,'position',s2pos);

    
    surf_chla=surf_chla(drange);
    total_diatC=total_diatC(drange);
    total_spC=total_spC(drange);
    total_zooC=total_zooC(drange);

    %[pathstr, fname, fext] = fileparts(filename); 
    %save_fig_figeps(fig,[fname, '_fig_oceiceline_dn']);
    
 
end
