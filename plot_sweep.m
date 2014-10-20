% Matlab function to plot sweeps from deevolve.
% J. Bloomfield (MIT) & J. Pearson (Durham)
% March 2014

% USEAGE: plot_sweep('LOGSDIR','RUNID');
% NOTE: 
% 1) "LOGSDIR" does not need a "slash"
% 2) "RUNID" does not have a ".dat"

% EXAMPLE:
% plot_sweep('logs','run0001');


function plot_sweep(varargin)
    close all;
    
    % Which data likelihoods do you want to plot?
    plotPLANCK = true;
    plotSNE = true;
    plotCOMBINED = false; % combined likelihood (lines only)
    plotCOMBINEDFILLED = true; % combined likelihood (filled contours)
    
    % Whats the font size for the plot label?
    fontsize_label = 17;
    % What units do you want to measure the physical size of the plot in?
    plot_sizeunits = 'inches';
    % Width of the plot (units of "plot_sizeunits")
    plot_width = 4.5;
    % Height of the plot (units of "plot_sizeunits")
    plot_height = 3.5;
    % Where will the plots go?
    plotdir = 'plots/';
    
    
    fprintf('Plots being output to %s\n',plotdir);
  
    
    % sigma-values
    sigma1 = (1-0.6827);
    sigma2 = (1-0.9545);
    sigma3 = (1-0.9973);
    
    % Get the data file
    datadir=char(varargin(1));
    RunID=char(varargin(2));
    filename = strcat(datadir,'/',RunID,'.dat');
    % It expects tab-delimited data
    delimiterIn = '\t';
    A = importdata(filename,delimiterIn);
    
    % The parameter and data combination names are stored at the top
    % We have to be careful getting the first parameter name: the data file
    % has a "#" at the start of the line
    p1name=char(A.textdata(1));
    p1name(1,1)=' ';
    p1name=strtrim(p1name);
    p2name=char(A.textdata(2));
    % Dump these parameter names into an array
    paramnames{1}=p1name;
    paramnames{2}=p2name;
    nparams=2;
    % Create latex strings for parameter names
    for p=1:nparams
        if strcmp(paramnames(p), 'wnaught')
            rep = '$w_0$';
            paramnames(p) = regexprep(paramnames(p),'wnaught',rep);
        end;
        if strcmp(paramnames(p), 'wa')
            rep = '$w_a$';
            paramnames(p) = regexprep(paramnames(p),'wa',rep);
        end;
        if strcmp(paramnames(p), 'Omegamh2')
            rep = '$\\Omega_{\\rm m} h^2$';
            paramnames(p) = regexprep(paramnames(p),'Omegamh2',rep);
        end;
         if strcmp(paramnames(p), 'Omegabh2')
            rep = '$\\Omega_{\\rm b} h^2$';
            paramnames(p) = regexprep(paramnames(p),'Omegabh2',rep);
        end;
        if strcmp(paramnames(p), 'Omegakh2')
            rep = '$\\Omega_{\\rm k} h^2$';
            paramnames(p) = regexprep(paramnames(p),'Omegakh2',rep);
        end;
    end;
    
    % Get arrays of the parameter values
    p1val=A.data(:,1);
    p2val=A.data(:,2);    
    
    % Get min & max of parameter values
    minp1=min(p1val);
    maxp1=max(p1val);
    minp2=min(p2val);
    maxp2=max(p2val);
    
    % Get the number of samples of each parameters
    % (used to construct the surface grid)
    x=p1val(1); x0=x;
    y=p2val(1);
    np1=0; np2=0;
    for i=2:size(p1val,1)
        if x==x0
            np1=0;
        end;
        if x~=p1val(i)
            x=p1val(i);
            np1=np1+1;
        end;
        if y~=p2val(i)
            y=p2val(i);
            np2=np2+1;
        end;
    end;
    
    % Put parameter values onto mesh for plotting
    [pp1,pp2]=meshgrid(linspace(minp1,maxp1,np1),linspace(minp2,maxp2,np2));
    
    % Get the different data-sets likelihoods
    WMAP=A.data(:,3);
    Planck=A.data(:,4);
    SN=A.data(:,5);
    combined=A.data(:,9);
    
    p1p2=[p1val,p2val];
    WMAP_D=griddata(p1val,p2val,WMAP,pp1,pp2);
    Planck_D=griddata(p1val,p2val,Planck,pp1,pp2);
    SN_D=griddata(p1val,p2val,SN,pp1,pp2);
  	combined_D=griddata(p1val,p2val,combined,pp1,pp2);
    
    hold on;
    if plotPLANCK
        contour(pp1,pp2,Planck_D,[sigma1],'LineWidth',3,'LineColor','b','DisplayName', 'Planck');
        contour(pp1,pp2,Planck_D,[sigma2],'LineWidth',2,'LineColor','b','HandleVisibility','off');
        contour(pp1,pp2,Planck_D,[sigma3],'LineWidth',1,'LineColor','b','HandleVisibility','off');
    end;
    if plotSNE
        contour(pp1,pp2,SN_D,[sigma1],'LineWidth',3,'LineColor','r','DisplayName', 'SNe');
        contour(pp1,pp2,SN_D,[sigma2],'LineWidth',2,'LineColor','r','HandleVisibility','off'); 
        contour(pp1,pp2,SN_D,[sigma3],'LineWidth',1,'LineColor','r','HandleVisibility','off');
    end;
    if plotCOMBINED
        contour(pp1,pp2,combined_D,[sigma1],'LineWidth',3,'LineColor','k','DisplayName', 'Combined');
        contour(pp1,pp2,combined_D,[sigma2],'LineWidth',2,'LineColor','k','HandleVisibility','off');
        contour(pp1,pp2,combined_D,[sigma3],'LineWidth',1,'LineColor','k','HandleVisibility','off');
    end;
    if plotCOMBINEDFILLED
        contourf(pp1,pp2,combined_D,[sigma1,sigma2,sigma3],'DisplayName', 'Combined');
    end;
    
    box on;
    legend show;
    legend(gca, 'boxoff');
    
    % Set limits on the plots
    xlim([minp1 maxp1]);
    ylim([minp2 maxp2]);
    
    % View contour plot from above
    view(0,90);
    % Put axes labels
    xlabel(paramnames(1),'interpreter','latex','fontsize',fontsize_label);
    ylabel(paramnames(2),'interpreter','latex','fontsize',fontsize_label);
     
    clear A;
    
    % Print plot to file
    outfile = strcat(plotdir,'plot_',RunID,'.pdf');
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',outfile);

    