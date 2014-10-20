function proc_1()

    close all;
    % directory for plots
    plotsDir = 'plots/';
    % file-name of the output plot
    ev_plots_fig = 'evPlots.pdf'; 
    Ct_plots_fig = 'Cplot.pdf';
    plot_sizeunits='inches';
    plot_width=8.5;
    plot_height=6.5;
    
    %The data IDs
    ID = 1;
    runIDs{ID} = '0001'; ID=ID+1;
    runIDs{ID} = '0002'; ID=ID+1;
    runIDs{ID} = '0004'; ID=ID+1;
    runIDs{ID} = '0006'; ID=ID+1;    
    runIDs{ID} = '0009'; ID=ID+1;    
    
    % legend labels
    v1 = '1';
    v2 = '1.1';
    v3 = '1.3';
    v4 = '1.5';
    v5 = '1.8';
    
    dataDir = 'logs/'; % where are all the files?
    runPrototype = 'run'; % what're the data files prefixed with?
    runSuffix = '.dat'; % what're the data files suffixed with?
    ev_plots_fig = strcat(plotsDir,ev_plots_fig);
    Ct_plots_fig = strcat(plotsDir,Ct_plots_fig);
    % column numbers in data file...
    tloc = 1;
    aloc = 2;
    zloc = 3;
    Hloc = 4;
    Hdotloc = 5;
    philoc = 6;
    phidotloc = 7;
    Omloc = 9;
    Orloc = 10;
    Odeloc = 12;
    wtotloc = 13;
    
    nIDs = size(runIDs);
    nIDs = nIDs(:,2); % number of files
    npx=4;  
    npy=2;
    xmx=0;
    cols={'g','r','b','k','m'};
    hold on;
    for n = 1:nIDs
        runID = runIDs{n};
        runFile = strcat(dataDir,runPrototype,runID,runSuffix);
        runData = load(runFile);

        % Get data    
        t = runData(:,tloc);
        a = runData(:,aloc);
        z = runData(:,zloc);
        H = runData(:,Hloc);
        Hdot = runData(:,Hdotloc);
        phi = runData(:,philoc);
        phidot = runData(:,phidotloc);
        Om = runData(:,Omloc);
        Or = runData(:,Orloc);
        Ode = runData(:,Odeloc);
        wde = runData(:,wtotloc);
        
        % construct Ricci scalar & integrate it
        nvals = size(a);
        nvals_use = nvals(1)-100;
        vol_int = 0.0;
        R_int = 0.0;
        smallest_C = 1E12;
        a_at_sC = 0.0;
        for i=1:nvals_use
            % The Ricci scalar
            RicciScalar = 6. * ( Hdot(i) + H(i) * H(i) ) / ( a(i) * a(i) );
            % Compute |\sqrt{-g}R|, and log-it
            logRicci(i) = log10(abs(RicciScalar));
            t_use(i) = t(i);
            volfac = a(i)*a(i)*a(i)*a(i);
            vol_int = vol_int + volfac;
            R_int = R_int + volfac*RicciScalar;
            int_integrand(i) = R_int / vol_int;
            if abs(int_integrand(i)) < smallest_C
                smallest_C = abs(int_integrand(i));
                a_at_sC = a(i);
                t_at_sC = t(i);
            end;
        end;
        Constraint_check = R_int/vol_int;
        fprintf('ID: %s, C = %f, a_last = %f, R_last = %f, ',runIDs{n},Constraint_check,a(nvals_use), logRicci(nvals_use));
        fprintf('C_min = %f, a(C_min) = %f, t(C_min) = %f\n',smallest_C,a_at_sC,t_at_sC);
       
        
        % get maximum value of t to use as xlim on plots    
        mint = min(t);
        maxt = t(nvals_use);
        pnum=1;
        if maxt>xmx
            xmx = maxt+0.5;
        end;
        
        % Do the plotting
        
        % 0 - C(\tau) :: value of <R> as a function of \tau_max
        hold on;
        %subplot(npx,npy,pnum);pnum=pnum+1;
        C_plot = figure(1);
        plot(t_use,int_integrand,'color',cols{n});
        refline(0,0);
        xlim([0 xmx]);
        ylim([-2E1 1E2]);
        xlabel('$\tau_{\rm max}$','interpreter','latex');
        ylabel('$\langle R\rangle(\tau_{\rm max}{\,})$','interpreter','latex');
        
        % 1 - scale factor
        hold on;
        ev_plots = figure(2);
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,a,'color',cols{n});
        refline(0,1);
        xlim([0 xmx]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$a(\tau)$','interpreter','latex');

        % 2 - Hubble
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,H,'color',cols{n});
        xlim([0 xmx]);
        ylim([-10 10]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$H$','interpreter','latex');
        
        % 3 - log10RicciMeasure
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t_use,logRicci,'color',cols{n});
        xlim([0 xmx]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$\log\,_{10}{}\, \left(\left|\sqrt{-g}{\,}R\right|\right)$','interpreter','latex');
        
        % 4 - H-dot
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,Hdot,'color',cols{n});
        xlim([0 xmx]);
        ylim([-10 2]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$\dot{H}$', 'interpreter','latex');
        
        % 5 - w_de
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,wde,'color',cols{n});
        box on;
        refline(0,1); % reference line: w = +1
        refline(0,-1); % reference line: w = -1
        xlim([0 xmx]);
        ylim([-1.1 2]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$w_{\rm tot}$', 'interpreter','latex');
        
        % 6 - OmegaDE
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,Ode,'color',cols{n});
        box on;
        xlim([0 xmx]);
        ylim([-0.1 1.1]);
        xlabel('$\tau$','interpreter','latex');
        ylabel('$\Omega_{\rm de}$', 'interpreter','latex');
        
        % 7 - phi
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,phi,'color',cols{n});
        xlim([0 xmx]);
        ylim([-5 2.1]);
        box on;
        xlabel('$\tau$','interpreter','latex');
        ylabel('$\phi$', 'interpreter','latex');
        
        % 8 - phidot
        hold on;
        subplot(npx,npy,pnum);pnum=pnum+1;
        plot(t,phidot,'color',cols{n});
        xlim([0 xmx]);
        ylim([ -1 0.1]);
        box on;
        xlabel('$\tau$','interpreter','latex');
        ylabel('$\dot{\phi}$', 'interpreter','latex');
        
        % Deallocate the Ricci & t_use arrays 
        %(else retains the same size which gives silly numbers)
        clear logRicci;
        clear t_use;
        clear int_integrand;
    end;
    % Setup legend
    L = legend(v1,v2,v3,v4,v5);
    set(L,'location','southwest');
    legend boxoff;
    hold off;
    % Plot to pdf
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',ev_plots_fig);
    
    C_plot_ax = findobj(C_plot,'type','axes');
    box(C_plot_ax);
    L = legend(C_plot_ax,v1,v2,v3,v4,v5);
    legend(C_plot_ax,'boxoff');
    
    
    set(C_plot, 'PaperUnits',plot_sizeunits);
    set(C_plot, 'PaperSize',[4 3]);
    set(C_plot, 'PaperPosition',[ 0 0 4 3]);
    set(C_plot, 'renderer', 'painters');
    print(C_plot,'-dpdf',Ct_plots_fig);
    
    
% eof    