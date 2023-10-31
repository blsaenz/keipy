% plt_nc -- modifed into function, Ben Saenz 7/2011

% !!! hardwired for hourly time step !!!!

% quick ice-ocean plots

function kei_plot_eco(kei,d)

	% hardwired for hourly time step
	time2=zeros(length(kei.step),1);
	dtday=1/24.; 
	stday=mod(kei.day(1),365);
	time2 = double(stday) + double(kei.step).*dtday;
	
	xlen=length(kei.zml);
	x=[0:xlen-1];
	
	figure;
	%pause;
	
	clf
	subplot(5,3,1)
	pcolor(time2,kei.zgrid(1:d),kei.T(1:d,:)); hold on;
    shading flat;
	plot(time2,-kei.hmx,'k','Linewidth',[2])
	plot(time2,kei.zml,'w')
	xlabel('Time (day)')
	ylabel('Depth')
	title('Temperature')
	colorbar
	
	subplot(5,3,4)
	pcolor(time2,kei.zgrid(1:d),kei.S(1:d,:)); hold on;
    shading flat;
	plot(time2,-kei.hmx,'k','Linewidth',[2])
	plot(time2,kei.zml,'w')
	xlabel('Time (day)')
	ylabel('Depth')
	title('Salinity')
	colorbar

    subplot(5,3,7)
	plot(time2,kei.hsn)
	ylabel('Thickness (m)')
	title('Snow Thickness')
    axx=axis;
    axx(1:2) = [min(time2) max(time2) ];
    axis(axx);
	subplot(5,3,10)
	plot(time2,kei.hi)
	ylabel('Thickness (m)')
	title('Ice Thickness')
    axx=axis;
    axx(1:2) = [min(time2) max(time2) ];
    axis(axx);
	subplot(5,3,13)
	plot(time2,kei.fice)
	xlabel('Time (day)')
	ylabel('Ice Fraction')
	title('Ice Fraction')
    axx=axis;
    axx(1:2) = [min(time2) max(time2) ];
    axis(axx);
	
    subplot(5,3,2)
	pcolor(time2,kei.zgrid(1:d),kei.NO3(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('NO3')
	colorbar

    subplot(5,3,5)
	pcolor(time2,kei.zgrid(1:d),kei.PO4(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('PO4')
	colorbar

    
    subplot(5,3,8)
	pcolor(time2,kei.zgrid(1:d),kei.SiO3(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('SiO3')
	colorbar

    subplot(5,3,11)
	pcolor(time2,kei.zgrid(1:d),kei.Fe(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('Fe')
	colorbar

    subplot(5,3,14)
	pcolor(time2,kei.zgrid(1:d),kei.NH4(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('NH4')
	colorbar

    subplot(5,3,3)
	pcolor(time2,kei.zgrid(1:d),kei.diatC(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('diatC')
	colorbar

    subplot(5,3,6)
	pcolor(time2,kei.zgrid(1:d),kei.diatChl(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('diatChl')
	colorbar

    
    subplot(5,3,9)
	pcolor(time2,kei.zgrid(1:d),kei.spC(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('spC')
	colorbar

    subplot(5,3,12)
	pcolor(time2,kei.zgrid(1:d),kei.zooC(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('zooC')
	colorbar

    subplot(5,3,15)
	pcolor(time2,kei.zgrid(1:d),kei.DIC(1:d,:)); hold on;
    shading flat;
	xlabel('Time (day)')
	ylabel('Depth')
	title('DIC')
	colorbar
    
    
    
end
