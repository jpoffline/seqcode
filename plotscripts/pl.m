function pl

    close all;

    dataDIR = '../testcut/';
    dataIDS = 'test22_2';
    dataSUF = '.dat';
    dataFIL = strcat(dataDIR,dataIDS,dataSUF);
    
    outfigname = dataIDS;
    
    data = load(dataFIL);
    size(data);
    
    plot_sizeunits='inches';
    plot_width=5;
    plot_height=3.5;
    axis_fontsize = 20;
    tick_fontsize = 15;
    cbaxis_fontsize = 15;
    
    p1loc = 1;
    p2loc = 2;
    HIRloc = 4;
    amaxloc = 7;
    w0loc = 8;
    tmaxloc = 9;
    likelog = 13;
    
    p1 = 10.^(data(:,p1loc));
    p2 = data(:,p2loc);
    HIR = log10(abs(data(:,HIRloc)));
    amax =  ( (data(:,amaxloc)));
    w0 =  10.^2*(1+( (data(:,w0loc))));
    tmax =   (data(:,tmaxloc));
    likelihood = exp(-0.5*data(:,likelog)/100);
    HIRv = data(:,HIRloc);
    
    np_x = 2;
    np_y = 2;
    np = 1;
    
    ref_m = - 1.1/10.^3;
    ref_c = - 1.8 / 10.^3;
    
    ss=size(HIR);
    
    h_thresh = 2.5;
    
    for i=1:ss(1)
        l = likelihood(i);
        h = HIR(i);
        if h > h_thresh
            l = 0.0;
        end;
        val(i) = -1;
        if HIRv(i) < 0
            val(i) = 1;
        end;
        
        likelihood(i) = l;
        
        p(i) = ref_m * p1(i) + ref_c;
        
    end;
    %size(HIRv)
    plot(p,likelihood);
    
    
    
   
    xmin=min(p1);
    xmax=max(p1);
    ymin=min(p2);
    ymax=max(p2);
    
    %{
    
    xlin = linspace(xmin, xmax, 1000);
    ylin = linspace(ymin, ymax, 1000);
    [X,Y] = meshgrid(xlin,ylin);

    %subplot(np_x,np_y,np); np = np + 1;
    Z = griddata(p1,p2,HIR,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    ylabel(cb, '$\log_{10}|\langle R\rangle|$','interpreter','latex','fontsize',cbaxis_fontsize);
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    ref = refline(ref_m,ref_c);
    set(ref,'LineStyle','--');
    set(ref,'LineWidth',5);
    set(ref,'color','black');
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_HIR.pdf'));
    
    %subplot(np_x,np_y,np); np = np + 1;
    Z = griddata(p1,p2,likelihood,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    ref = refline(ref_m,ref_c);
    set(ref,'LineStyle','--');
    set(ref,'LineWidth',5);
    set(ref,'color','black');
    
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    ylabel(cb, '$\mbox{Likelihood}$','interpreter','latex','fontsize',cbaxis_fontsize);
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_L.pdf'));
    
   % subplot(np_x,np_y,np); np = np + 1;
    Z = griddata(p1,p2,amax,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    ylabel(cb, '$a_{\rm max}$','interpreter','latex','fontsize',cbaxis_fontsize);
    ref = refline(ref_m,ref_c);
    set(ref,'LineStyle','--');
    set(ref,'LineWidth',5);
    set(ref,'color','black');
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_amax.pdf'));
    
   % subplot(np_x,np_y,np); np = np + 1;
    Z = griddata(p1,p2,tmax,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    ylabel(cb, '$t_{\rm max}$','interpreter','latex','fontsize',cbaxis_fontsize);
    ref = refline(ref_m,ref_c);
    set(ref,'LineStyle','--');
    set(ref,'LineWidth',5);
    set(ref,'color','black');
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_tmax.pdf'));
    
    
    
    
    Z = griddata(p1,p2,w0,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    ylabel(cb, '$10^2(1+w_0)$','interpreter','latex','fontsize',cbaxis_fontsize);
    ref = refline(ref_m,ref_c);
    set(ref,'LineStyle','--');
    set(ref,'LineWidth',5);
    set(ref,'color','black');
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_w0.pdf'));
    
    
    
    Z = griddata(p1,p2,val,X,Y,'cubic');
    imagesc(xlin, ylin,Z);
    cb = colorbar;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    xlabel('$m^3$','interpreter','latex','fontsize',axis_fontsize);
    ylabel('$\Omega_{\rm k} h^2$','interpreter','latex','fontsize',axis_fontsize);
    ylabel(cb, '$\mbox{sign}\langle R\rangle$','interpreter','latex','fontsize',cbaxis_fontsize);
    set(gca,'YDir','normal');
    set(gca,'FontSize',tick_fontsize)
    
    set(gcf, 'PaperUnits',plot_sizeunits);
    set(gcf, 'PaperSize',[plot_width plot_height]);
    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
    set(gcf, 'renderer', 'painters');
    print('-dpdf',strcat(outfigname,'_HIRv.pdf'));
    
    %}
    