function like_1d

    close all;
    
    dataDIR = '../testcut/';
    dataIDS = 'test22_33';
    dataSUF = '.dat';
    dataFIL = strcat(dataDIR,dataIDS,dataSUF);
    
    data = load(dataFIL);
    
    
    plot_sizeunits='inches';
    plot_width=5;
    plot_height=3.5;
    axis_fontsize = 20;
    tick_fontsize = 15;
    cbaxis_fontsize = 15;
    outfigname = dataIDS;
    m = 10.^(data(:,1));
    Omk = data(:,2);
    L = data(:,4);
    amax = data(:,5);
    tmax = data(:,6);
    L = L / max(L);
    
    m_max = max(m);
    m_min = min(m);
    m_max = 0.75;
    m_min = 0.55;
    npx = 2;
    npy = npx;
    np = 1;
    s1 = 0.6;
    s2 = 0.9;
    
   % subplot(npx,npy,np); np = np+1;
    plot(m, L);
     xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\mbox{Likelihood}$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'FontSize',tick_fontsize);
    xlim([m_min m_max]);
    
     set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_L.pdf'));
    
    close all;
    
    %subplot(npx,npy,np); np = np+1;
     for ii=1:length(m)-1
        c = 'g';
        if(L(ii)) > s1
            c = 'b';
        end;
        if(L(ii)) > s2
            c = 'r';
        end;
        
        plot([m(ii) m(ii+1)],[amax(ii) amax(ii+1)],c,'linewidth',4);
        hold on;
    end;
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$a_{\rm max}$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'FontSize',tick_fontsize);
    xlim([m_min m_max]);
    
     set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_amax.pdf'));
    
    close all;
    %subplot(npx,npy,np); np = np+1;
    for ii=1:length(m)-1
        c = 'g';
        if(L(ii)) > s1
            c = 'b';
        end;
        if(L(ii)) > s2
            c = 'r';
        end;
        
        plot([m(ii) m(ii+1)],[tmax(ii) tmax(ii+1)],c,'linewidth',4);
        hold on;
    end;
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$t_{\rm max}$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'FontSize',tick_fontsize);
    xlim([m_min m_max]);
    
     set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_tmax.pdf'));
    
    close all;
   % subplot(npx,npy,np); np = np+1;
     for ii=1:length(m)-1
        c = 'g';
        if(L(ii)) > s1
            c = 'b';
        end;
        if(L(ii)) > s2
            c = 'r';
        end;
        
        plot([m(ii) m(ii+1)],[Omk(ii) Omk(ii+1)],c,'linewidth',4);
        hold on;
    end;
     xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k}h^2 (\mbox{derived})$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'FontSize',tick_fontsize);
    xlim([m_min m_max]);
     set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_Omk.pdf'));
    