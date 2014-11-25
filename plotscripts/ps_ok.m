function ps

    close all;
    axes_fontsize = 15;
    
    datadir = '../testcut/';
    datafil = 'test13.dat';
    
    data = load(strcat(datadir, datafil));
    
    omkloc = 2;
    HRloc = 4;
    amaxloc = 7;
    wnowloc = 8;
    tfinalloc = 9;
    
    
    omk = 1E3*data(:,omkloc);
    amax = data(:,amaxloc);
    wnow = log10((1+ data(:,wnowloc)));
    tfinal = data(:,tfinalloc);
    
    spacing = 0.03;
    left=0.1;
    bottom1=0.7;
    bottom2=0.4;
    bottom3=0.1;
    width=0.8;
    height=bottom1-bottom2-spacing; % which is also bottom1-bottom2
    
    
    
    axes('Position',[left bottom1 width height]);
    plot(omk, tfinal);
    hy1=ylabel('$t_{\rm final}$','interpreter','latex');
    set(gca, 'XTickLabel', []);
    
    
    axes('Position',[left bottom2 width height]);
    plot(omk, amax);
    hy2=ylabel('$a_{\rm max}$','interpreter','latex');
    set(gca, 'XTickLabel', []);
    
    axes('Position',[left bottom3 width height]);
    plot(omk, wnow);
    hx3=xlabel('$10^3\Omega_kh^2$','interpreter','latex');
    hy3=ylabel('$\log{}_{10}{\,}(1 + w)$','interpreter','latex');
    
    
    set(hx3, 'fontsize', axes_fontsize);
    set(hy3, 'fontsize', axes_fontsize);
    set(hy2, 'fontsize', axes_fontsize);
    set(hy1, 'fontsize', axes_fontsize);
    
     outfile = 'p2.pdf';
     % What units do you want to measure the physical size of the plot in?
    plot_sizeunits = 'inches';
    % Width of the plot (units of "plot_sizeunits")
    plot_width = 6.5;
    % Height of the plot (units of "plot_sizeunits")
    plot_height = 4.5;
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',outfile);
    