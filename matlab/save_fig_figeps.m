function save_fig_figeps(ghandle,fname)

   set(gcf, 'PaperPositionMode', 'auto');
   fname_ext = strcat(fname,'.fig');
   saveas(ghandle,fname_ext);
   fname_ext = strcat(fname,'.eps');
   %saveas(ghandle,fname_ext,'epsc2');
   print(ghandle,'-painters','-depsc',fname_ext);
end