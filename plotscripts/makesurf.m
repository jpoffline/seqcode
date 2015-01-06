function makesurf(varargin)

    close all;
    
    name = char(varargin(1));
    which = char(varargin(2));
    whichID = char(varargin(3));
    
    plot_sizeunits='inches';
    plot_width=4.5;
    plot_height=3.5;
   
    
   
    % where is the root of all output files?
    rootdir='../output/';
    % where will the plots go?
    plotdir = '../plots';
    % check plots directory exists, and if not: make it
    if ~exist(plotdir,'dir')
        mkdir(plotdir);
    end;
    % add on a slash to the plotdir name
    plotdir=strcat(plotdir,'/');
    
   
    infilename = strcat(rootdir,name,'/file_',which,'.dat');
    inputdata = load(infilename);
    outfigname = strcat(plotdir,name,'_',whichID,'_',which,'.pdf'); 
    
    x = inputdata(:,1);
    y = inputdata(:,2);
    xmin=min(x);
    xmax=max(x);
    ymin=xmin;
    ymax=xmax;
    axlim=20;
    scale=1.0;
    if strcmp(whichID,'cham')
        ID = 3;
    end;
    if strcmp(whichID,'chamforce')
        ID = 6;
    end;
    if strcmp(whichID,'gravforce')
        ID = 7;
        scale=1000;
    end;
    if strcmp(whichID,'phierr')
        ID = 8;
    end;
    datF = inputdata(:,ID)/scale;
    
    xlin = linspace(xmin, xmax, 1000);
    ylin = linspace(ymin, ymax, 1000);
    [X,Y] = meshgrid(xlin,ylin);

    Z = griddata(x,y,datF,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    ylabel(cb, whichID);
    colormap Hot;
    shading interp;
    xlabel('x');
    ylabel('y');
    
    xlim([-axlim axlim]);
    ylim([-axlim axlim]);
    caxis([0.5*max(datF) max(datF)])
    box on;
    
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',outfigname);