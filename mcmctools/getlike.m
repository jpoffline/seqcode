% Matlab function to plot MCMC output from deevolve.
% J. Bloomfield (MIT) & J. Pearson (Durham)
% March 2014
%
% USEAGE: getlike('run')
%
% This works by reading in the chain data, binning it,
% and smoothing the bins. 
% We then plot all 1D and 2D likelihoods - 1D go along the diagonal. 
	
function getlike(varargin)

    close all;
    % Where will the plots be stored?
    plotdir = 'plots/';
    % Whats the name of the priors file?
    PriorsFileName = 'priors.txt'; 
    % Should we plot the surface?
    plotsurf = true;
    % Should we plot contours?
    plotcontours = true;
    % Should we plot a random selection of the samples?
    plotsamples = false;
    % How many of the random samples should we plot?
    nrandsamples = 500;
    % How many bins?
    nbins = 45;
    
    contourlinewidth = 3;
    % Smoothing parameter for 1D likelihood
    smoothparam_1D = 5;
    % Smoothing parameter for 2D likelihood
    smoothparam_2D = 4;
    % Set the maximum number of iterations the smoother will perform
    smoothiternumbs = 500;
	
	% When making the plot, what units to use?
    plot_sizeunits='inches';
	% How wide is the plot?
	plot_width=15;
	% How high is the plot?
	plot_height=15;
	
    sigma1 = (1-0.6827);
    sigma2 = (1-0.9545);
    sigma3 = (1-0.9973);
    
    RunDirID = char(varargin);
    RunDir = strcat('../chains/',RunDirID,'/');
    RunDirContents = dir(fullfile(RunDir,'*chain*'));
    dims = size(RunDirContents);
    nchains = dims(1);
   
    priorsinfo = fopen(strcat(RunDir,PriorsFileName));
    C = textscan(priorsinfo, '%s %s %f %f %f', 'CommentStyle', '#');
    fclose(priorsinfo);
    paramsectn = C{:,1};
    paramnames = C{:,2};
    paramlower = C{:,3};
    paramupper = C{:,4};
    paramsigma = C{:,5};
    nparams = numel(paramnames)
    paramnames
    % Setup formatting stuff
    formatspec = '%f';
	% This is where labels from the code and priors file
	% can be converted into LaTeX labels.
    for p=1:nparams
        formatspec = strcat('%f',formatspec);
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
        if strcmp(paramnames(p), 'desiredh')
            rep = '$h$';
            paramnames(p) = regexprep(paramnames(p),'desiredh',rep);
        end;
        if strcmp(paramnames(p), 'phi0')
            rep = '$\\phi_0$';
            paramnames(p) = regexprep(paramnames(p),'phi0',rep);
        end;
        if strcmp(paramnames(p), 'phidot0')
            rep = '$\\dot{\\phi}_0$';
            paramnames(p) = regexprep(paramnames(p),'phidot0',rep);
        end;
        if strcmp(paramnames(p), 'mass')
            rep = '$m^3_{\\rm slope}$';
            paramnames(p) = regexprep(paramnames(p),'mass',rep);
        end;
    end;
   
    
   % Now open up the chains to compute bins etc
    for n = 1:nchains
        
        % Feed back some progress info to screen
        fprintf('Chain # %i of %i\n',n, nchains);
        
        % Open up this chain
        chainfilename = strcat(RunDir,RunDirContents(n).name);
        chainfile = fopen(chainfilename);
        chaindata = textscan(chainfile,formatspec,'CommentStyle', '#');
        fclose(chainfile);
        chaindata = cell2mat(chaindata);
        
        % Get some important dimensions info
        chaindims = size(chaindata);
        nsamples = chaindims(1)-1;
        %nparams = chaindims(2)-1;
        chainID = 1E4+n-1;
        
        plotnumber = 1;
        skip = 0;
        % Get the samples of parameters
        for p1=1:nparams
            
            % Uncomment this to find out which parameter is being looked at
            %fprintf('%s\n',char(paramnames(p1)))
            
            % Get the samples data
            param1sample=chaindata(:,p1);
            % Construct a histogram of this parameter
            [bindata1,w] = hist(param1sample,nbins);
            % Smoothe the bins
            smdat1 = smoothn(bindata1,smoothparam_1D);
            % Plot
            subplot(nparams,nparams,plotnumber);
            plot(w,smdat1);
			ylim([min(smdat1) max(smdat1)]);
            set(gca, 'YTick', []);
            xlabel(char(paramnames(p1)), 'Interpreter','Latex');
            
            plotnumber=plotnumber+1;
            
            for p2=p1+1:nparams
                  
                 % Get the samples for this parameter
                 param2sample=chaindata(:,p2);   
                 
                 % Get the indexing of the subplot for this 2D likelihood plot
                 subplot(nparams,nparams,plotnumber);
                 
                 % Construct a 2D data array
                 DATA = [param1sample, param2sample];
                 
                 % Obtain smoothed histogram of the 2D data
                 [H,X,Y] = smoothhist2D(DATA, smoothparam_2D, [nbins nbins],0.1*sigma3);
                 colormap jet;
                 
                 if plotsurf
                    hold on;
                 end;
                 
                 % Plot contours
                 if plotcontours
                    maxbin = max(max(H));  
                    contour(X,Y,H,sigma1*maxbin,'LineWidth',2,'Linecolor','k');
                    hold on;
                    contour(X,Y,H,sigma2*maxbin,'LineWidth',1,'Linecolor','k');
                    hold on;
                 end;
                 
                 %  Plot (randomly) from the samples
                 if plotsamples
                     freezeColors;
                     % First, shuffle the rows randomly;
                     shuffsamps = randperm(nsamples);
                     % Create array from the "top" howmanys of the shuffled
                     % samples
                     randsamps = chaindata(shuffsamps(1:nrandsamples),:);
                     % Create array holding all of the randomly sampled
                     % parameters of type "p1"
                     Xrs = randsamps(:,p1);
                     % and of type "p2"
                     Yrs = randsamps(:,p2);
                     % Get the normalised likelihoods
                     randlike = randsamps(:,nparams+1) / sum(randsamps(:,nparams+1));
                     % Set the size of the points
                     S = 1;
                     % Get the color of the point based on the value of the
                     % likelihood
                     C = randlike;
                     % Plot in scatter diagram
                     
                     scatter(Xrs,Yrs,S,C,'fill');
                     colormap autumn;
                     % Flip the colormap
                     colormap(flipud(colormap))
                     freezeColors;
                     hold on;
                 end;
                 
                 
                 % Set the limits on the plot; 
				 % inherit from the smoothed histogram
                 xlim([min(X) max(X)]);
                 ylim([min(Y) max(Y)]);
                 
                 % Put the x & y labels on
                 xlabel(char(paramnames(p1)), 'Interpreter','Latex');
                 ylabel(char(paramnames(p2)), 'Interpreter','Latex');
                 
                 hold off;
                 % Increment plotnumber counter
                 plotnumber = plotnumber + 1;
                 
            end;
            skip = skip + 1;
            plotnumber = plotnumber + skip;
            
        end; % finish looking at all parameters
        
        % Create plot & print as pdf
        outfile = strcat(plotdir,'tri_',RunDirID,'_',num2str(chainID),'.pdf');
        set(gcf, 'PaperUnits',plot_sizeunits);
	    set(gcf, 'PaperSize',[plot_width plot_height]);
	    set(gcf, 'PaperPosition',[ 0 0 plot_width plot_height]);
	    set(gcf, 'renderer', 'painters');
	    print('-dpdf',outfile);
        fprintf('Plot saved %s\n',outfile);
    end; % finish reading in all the chains
    
    % Print a final message
    
    fprintf('Done\n');
    