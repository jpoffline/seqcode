function plotdistsurf(Border,Sigma,stepSize,Data)

    X = Data(:,1);
    Y = Data(:,2);
    dims = size(Data);
    
    NumData = length(X);
    

    Xrange = [min(X)-Border max(X)+Border];
    Yrange = [min(Y)-Border max(Y)+Border];


    %Setup coordinate grid
    [XX YY] = meshgrid(Xrange(1):stepSize:Xrange(2), Yrange(1):stepSize:Yrange(2));
    YY = flipud(YY);

    %Parzen parameters and function handle
    pf1 = @(C1,C2) (1/NumData)*(1/((2*pi)*Sigma^2)).*...
         exp(-( (C1(1)-C2(1))^2+ (C1(2)-C2(2))^2)/(2*Sigma^2));

    PPDF1 = zeros(size(XX));    

    %Populate coordinate surface
    [R C] = size(PPDF1);
    NN = length(Data);
    for c=1:C
    for r=1:R 
       for d=1:NumData 
            PPDF1(r,c) = PPDF1(r,c) + ...
                pf1([XX(1,c) YY(r,1)],[Data(d,1) Data(d,2)]); 
       end
    end
    end


    %Normalize data
    m1 = max(PPDF1(:));
    PPDF1 = PPDF1 / m1;

    %Set up visualization
    set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',20)
    fig = figure(1);clf
    stem3(Data(:,1),Data(:,2),zeros(NumData,1),'b.');
    hold on;

    %Add PDF estimates to figure
    s1 = surfc(XX,YY,PPDF1);shading interp;alpha(s1,'color');
    sub1=gca;
    view(2)
    axis([Xrange(1) Xrange(2) Yrange(1) Yrange(2)])
