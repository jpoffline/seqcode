function scan()

    close all;

    
    
    dataDir = 'logs/';
    runPrototype = 'run';
    runDataSuffix = '.dat';
    runLogSuffix = '.log';
    ntrail = 4; % number of trailing zeros in file name
    
    logs = dir(strcat(dataDir,'*',runLogSuffix));
    ss = size(logs);
    nlogs = ss(1);
    
    
    param_min = 5.5;
    param_max = 6.0;
    dparam = 0.005;
    
    npx=2;
    npy=1;
    
    tloc = 1;
    aloc = 2;
    zloc = 3;
    Hloc = 4;
    Hdotloc = 5;
    
    hold on;
    for n = 1:nlogs
        np=1;
        param = param_min + n *dparam;
        % Get the file name
        filnametemp = int2str(n);
        while size(filnametemp) < ntrail
            filnametemp = strcat('0',filnametemp);
        end
        
        runID = filnametemp;
        runFile = strcat(dataDir,runPrototype,runID,runDataSuffix);
        runData = load(runFile);

        t = runData(:,tloc);
        a = runData(:,aloc);
        z = runData(:,zloc);
        H = runData(:,Hloc);
        Hdot = runData(:,Hdotloc);
        nvals = size(H);
        nvals_use = nvals(1)-4;

        vol_int = 0.0;
        R_int = 0.0;
        for i=1:nvals_use
            RicciScalar = 6. * ( Hdot(i) + H(i) * H(i) ) / ( a(i) * a(i) );
            logRicci(i) = ((RicciScalar));
            t_use(i) = t(i);
            volfac = a(i)*a(i)*a(i)*a(i);
            vol_int = vol_int + volfac;
            R_int = R_int + volfac*RicciScalar;
        end;
        Constraint_check = R_int/vol_int;
        fprintf('ID:%s, P = %f, C = %f\n', runID, param, Constraint_check);
        C(n) = Constraint_check;
        
        col = [255*n/nlogs/255,0,0];
        
        
        hold on;
        subplot(npx,npy,np); np = np+1;
        plot(t_use,logRicci,'color',col); 
        ylabel('log_{10} R(\tau)');
        
        xlim([7 7.2]);
       % ylim([0 1E3]);
        
        
        
        
    end;
    
    m = [param_min:dparam:param_max];
    
    subplot(npx,npy,np);
    plot(m,C);
    refline(0,0);
    xlabel('m^3_{slope}');
    ylabel('Historic average of Ricci scalar');