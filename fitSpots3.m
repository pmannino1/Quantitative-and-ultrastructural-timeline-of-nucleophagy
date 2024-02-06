function [new_spotData, params, FitResults_G, FitResults_GE] = fitSpots3(Image, spotData, varargin)

% very much in construction
% Fit spots in provided images given initial estimates of spots of
% interest

% check what input is provided as a spotList
% for now we will assume that output of findIrregularSpots3 is provided



%specify default settings in case the user specifies no additional values
fitRadius = 9;
marginRadius = 10;
imageScale=[0,1];
markerSize=4;
fitWeight='srqt';
showFitData=true;

version='v0.1: 2021.06.04';

for k = 1:length(varargin)
    
    if strcmpi(varargin{k},'fitRadius')
        fitRadius = varargin{k+1};
    elseif strcmpi(varargin{k},'fitWeight')
        fitWeight = varargin{k+1};
    elseif strcmpi(varargin{k},'marginRadius')
        marginRadius = varargin{k+1};
    elseif strcmpi(varargin{k},'imageScale')
        imageScale = varargin{k+1};
    elseif strcmpi(varargin{k},'markerSize')
        markerSize = varargin{k+1};
     elseif strcmpi(varargin{k},'showFitData')
        showFitData = varargin{k+1};
        %     elseif strcmpi(varargin{k},'quantileThreshold')
        %         quantileThreshold = varargin{k+1};
        %     elseif strcmpi(varargin{k},'intensityRatioThreshold')
        %         intensityRatioThreshold = varargin{k+1};
    end
end


% image cropping params
mrgn=marginRadius; % margin around the spot center  for the image cropping
% spot fitting params
R_fit=fitRadius;    % radius of the circular area arounfd the spot center to be used for fitting

% VISUALIZATION params
im_sc= imageScale; % contrast, in fractions of dynamic range of provided image (i.e. relAtive to max and min)
ms=markerSize; % marker size for the 3D image/fit plot

% For future to have different inputs workable...
% but currently only one
input_type='fIrS3';


% PREPARATION of auxilaury variables etc
xx=1:2*mrgn+1;
XX=repmat(xx, 2*mrgn+1,1);
YY=XX';
RR2=(XX(:)-mrgn-1).^2+(YY(:)-mrgn-1).^2;


% check weights options
    switch fitWeight
        case 'unity'
            disp('using "unity" weights for fitting')
       case 'sqrt'
           disp('using "sqrt" weights for fitting')
        case 'linear'
            disp('using "linear" weights for fitting');
        otherwise
            disp('Sorry such option is not provided')
            disp('currently only "sqrt"(default), "unity", "linear"  are available ')
            disp('using "unity" weights now')
    end

% generating fig windows
fig_spot3D=figure;
fig_spot1D_0=figure;
fig_spot1D_1=figure;
fig_spot=figure;

% prepare collectors for fit results
FitResults_G=[];
FitResults_GE=[];

%% SPOT FITTING  one-by-one from specified list

n_spots=length(spotData);

for ii=1:n_spots
    
    %frame=spots_2_fit(ii, 1);
    %spotID=spots_2_fit(ii, 2);
    
    % get spot data
    switch input_type
        case 'fIrS3'
            spotData_1=spotData{ii};
            spotPos=spotData_1.spotPosition;
            x=spotPos(1);
            y=spotPos(2);
            z=spotPos(3);
            pos=round([x,y]);
            frame=round(z);
        otherwise
            disp('Sorry, do not know what to do with data provided as a spotList')
            disp('spotData is expected to be an ouput from findIrregularSpots3')
            return
    end
    
    % CROP ROI from the image
    % pad image in case spots is close to the border
    image_z=Image(:,:,frame);
    im_size=size(image_z);
    im_class=class(image_z);
    imageX1=eval([im_class,'(zeros([',num2str([im_size(1)+mrgn*2,im_size(2)+mrgn*2]),']));']);
    im_sizeX=size(imageX1);
    imageX1(mrgn+1:im_sizeX(1)-mrgn, mrgn+1:im_sizeX(2)-mrgn)=image_z;
    
    % crop ROI
    image_1=imageX1(pos(2)+mrgn-mrgn:pos(2)+mrgn+mrgn, pos(1)-mrgn+mrgn:pos(1)+mrgn+mrgn);
    im_min=min(image_1(image_1>0));
    im_max=max(image_1(image_1>0));
    im_range=im_max-im_min;
    im_0=im_min+im_sc(1)*im_range;
    im_1=im_min+im_sc(2)*im_range;
    
    % set weights
    switch fitWeight
        case 'unity'
            WW= ones(size(XX));
       case 'sqrt'
            WW= sqrt(double(image_1(:)-min(image_1(:))));
        case 'linear'
            WW= double(image_1(:)-min(image_1(:)));
        otherwise
            WW= ones(size(XX));
    end
    
    % FITTING
    % define fit function:
    % Gaussian +constant
    model_G = fittype('c+ b*exp(-(x-x0).^2/(2*sigma2)-(y-y0).^2/(2*sigma2))',...
        'dependent',{'z'},'independent',{'x','y'},...
        'coefficients',{'c',  'x0', 'y0',  'b', 'sigma2'});
    % Gaussian+ erf for backround
    model_GE = fittype('c+a*erf((cos(phi)*(x-x0)-sin(phi)*(y-y0))/gamma) + b*exp(-(x-x0).^2/(2*sigma2)-(y-y0).^2/(2*sigma2))',...
        'dependent',{'z'},'independent',{'x','y'},...
        'coefficients',{'c', 'a','phi', 'x0', 'y0', 'gamma', 'b', 'sigma2'});
% model_GE = fittype('c+a*erf((cos(phi)*(x-x0)-sin(phi)*(y-y0))/gamma) + b*exp(-(x-x0).^2/(2*sigma2)-(y-y0).^2/(2*sigma2)) + a*()erf((cos(phi)*(x-x0)-sin(phi)*(y-y0))/gamma) + b*exp(-(x-x0).^2/(2*sigma2)-(y-y0).^2/(2*sigma2))',...
%         'dependent',{'z'},'independent',{'x','y'},...
%         'coefficients',{'c', 'a','phi', 'x0', 'y0', 'gamma', 'b', 'sigma2'});
        
    % intial estimates for coefficients:
    coeff_G_0= [min(image_1(:)), mrgn+1, mrgn+1, max(image_1(:))-median(image_1(:)), 3];
    %coeff_GE_0= [min(image_1(:)), std(double(image_1(:))), 1, mrgn+1, mrgn+1, 2, max(image_1(:))-median(image_1(:)), 3];
    %coeff_GE_0= [min(image_1(:)), 2*(max(image_1(:))-median(image_1(:))), 1, mrgn+1, mrgn+1, 2, max(image_1(:))-median(image_1(:)), 3];
    coeff_GE_0= [min(image_1(:)), std(double(image_1(:))), 1, mrgn+1, mrgn+1, 2, max(image_1(:))-median(image_1(:)), 3];
    
    % set bounadries for coefficeints:
    lb_G=[0, 0, 0, 0, 0];
    ub_G=[2*max(image_1(:)), 2*mrgn+1, 2*mrgn+1,  Inf, 12.25];
    lb_GE=[0, 0, 0, 0, 0, 0, 0, 0];
    ub_GE=[2*max(image_1(:)), Inf, 2*pi, 2*mrgn+1, 2*mrgn+1, 2*mrgn, Inf, 12.25];
    
    % convert px value into double as Matlab don't like integers sometime:
    ZZ=double(image_1(:));
    
    % set weights
    switch fitWeight
        case 'unity'
            WW= ones(size(XX));
       case 'sqrt'
            WW= sqrt(ZZ-min(ZZ));
        case 'linear'
            WW= ZZ-min(ZZ);
        otherwise
            disp('Sorry such option is not provided')
            disp('currently only "sqrt"(default), "unity", "linear"  are available ')
            disp('using "unity" weights now')
            WW= ones(size(XX));
    end
    
    % fit only pixel witin radius R_fit
    ind=RR2<R_fit^2;
    
    % Finally fitting
    [fit_0, gof_0, fit_out_0] = fit([XX(ind), YY(ind)], ZZ(ind),model_G,'Startpoint',coeff_G_0,'Lower',lb_G,'Upper',ub_G, 'Weights', WW(ind));
    [fit_1, gof_1, fit_out_1] = fit([XX(ind), YY(ind)], ZZ(ind),model_GE,'Startpoint',coeff_GE_0,'Lower',lb_GE,'Upper',ub_GE, 'Weights', WW(ind));
    if showFitData
        disp(fit_0)
        disp(fit_1)
    end
     
    % SHOW spot
    % show as an image
    figure(fig_spot);
    set(gcf,'Name',['frame = ',num2str(frame), ' spot # ',num2str(ii)])
    hold off
    imshow(image_1,[im_0,im_1],'InitialMagnification',1600,'Border','tight');
    hold on
    plot(mrgn+1+x-pos(1), mrgn+1+y-pos(2), 'xc')
    set(gcf,'Name',['frame = ',num2str(frame), ' spot # ',num2str(ii)])
    
    % show image and fit as 3D plot
    figure(fig_spot3D);
    hold off
    plot3(XX(:), YY(:), image_1(:), 'o', 'MarkerSize', ms, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
    hold on
    
    
    % Collect and show fit results
    % fit_0 - Gauss+C
    fit_X=fit_0;
    gof_X=gof_0;
    col_X=[0.5, 0, 0.5];
    % get coeff values, conf intervals and some goodness of fit measures
    coeff=coeffvalues(fit_X);
    c=coeff(1);
    x_0=coeff(2);
    y_0=coeff(3);
    b=coeff(4);
    sigma2=coeff(5);
    
    coeff_int=confint(fit_X);
    delta_coeff_int=diff(coeff_int,1);
    av_err_sum=sum(delta_coeff_int./coeff)/length(coeff);
    adjRsq=gof_X.adjrsquare;
    rmse=gof_X.rmse;
    
    FitResults_G.coeff(ii, :)=coeff;
    FitResults_G.conf_int(ii, :)=reshape(coeff_int, 1, 2*size(coeff_int,2));
    FitResults_G.GoF(ii, :)=[av_err_sum, adjRsq, rmse];
    FitResults_G.spotInfo(ii, :)=[frame, ii];
    
    % show best-fit spot center
    figure(fig_spot);
    plot(x_0, y_0, '+','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',col_X);    
    
    % show best-fit surface
    figure(fig_spot3D)
    ZZfit = c+ b*exp(-(XX-x_0).^2/(2*sigma2)-(YY-y_0).^2/(2*sigma2));
    mesh(XX, YY, ZZfit, 'EdgeColor', col_X);
    
    % show image and fit as 1D slice through the center of the fit
    figure(fig_spot1D_0);
    hold off
    ind_x= (XX==round(x_0));
    ind_y= (YY==round(y_0));
    % plot(xx, image_1(round(x_0),:),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',colX);
    plot(xx, image_1(ind_x),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',col_X);
    hold on
    plot(xx, image_1(ind_y),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor','w');
    plot(xx, ZZfit(ind_x), '-', 'LineWidth',1, 'Color', col_X);
    plot(xx, ZZfit(ind_y),':', 'LineWidth',1, 'Color', col_X);
    
    %fit_1 - Gauss+Erf
    fit_X=fit_1;
    gof_X=gof_1;
    col_X=[0.5, 0., 0];
    % get coeff values, conf intervals ans smoe goodness of fit measures
    coeff=coeffvalues(fit_X);
    c=coeff(1);
    a=coeff(2);
    phi=coeff(3);
    x_0=coeff(4);
    y_0=coeff(5);
    gamma=coeff(6);
    b=coeff(7);
    sigma2=coeff(8);
    
    coeff_int=confint(fit_X);
    delta_coeff_int=diff(coeff_int,1);
    av_err_sum=sum(delta_coeff_int./coeff)/length(coeff);
    adjRsq=gof_X.adjrsquare;
    rmse=gof_X.rmse;
    
    FitResults_GE.coeff(ii, :)=coeff;
    FitResults_GE.conf_int(ii, :)=reshape(coeff_int, 1, 2*size(coeff_int,2));
    FitResults_GE.GoF(ii, :)=[av_err_sum, adjRsq, rmse];
    FitResults_GE.spotInfo(ii, :)=[frame, ii];
    
     % show best-fit spot center
    figure(fig_spot);
    plot(x_0, y_0, '+','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',col_X);    
   
     % show best-fit surface
    figure(fig_spot3D)
    ZZfit = c+a*erf((cos(phi)*(XX-x_0)-sin(phi)*(YY-y_0))/gamma) + b*exp(-(XX-x_0).^2/(2*sigma2)-(YY-y_0).^2/(2*sigma2));
    mesh(XX, YY, ZZfit, 'EdgeColor', col_X);
    
        % show image and fit as 1D plot
    figure(fig_spot1D_1);
    hold off
    ind_x= (XX==round(x_0));
    ind_y= (YY==round(y_0));
    % plot(xx, image_1(round(x_0),:),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',colX);
    plot(xx, image_1(ind_x),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor',col_X);
    hold on
    plot(xx, image_1(ind_y),'o','MarkerSize',ms,'MarkerEdgeColor',col_X, 'MarkerFaceColor','w');
    plot(xx, ZZfit(ind_x), '-', 'LineWidth',1, 'Color', col_X);
    plot(xx, ZZfit(ind_y),':', 'LineWidth',1, 'Color', col_X);   
    
    figure(fig_spot3D);
     legend({'exp image', 'Gauss + c', 'Gauss + erf'}, 'Location', 'northeast' , 'FontSize', 10)
     set(gcf,'Name',['frame = ',num2str(frame), ' spot # ',num2str(ii)])
        
    figure(fig_spot1D_0);
     legend({'X-slice', 'Y-slice', 'Fit X-slice', 'Fit Y-slice'}, 'Location', 'northeast' , 'FontSize', 10)
     set(gcf,'Name',['frame = ',num2str(frame), ' spot # ',num2str(ii)])
     title ('Gauss + c fit', 'FontSize', 10)
 
    figure(fig_spot1D_1);
     legend({'X-slice', 'Y-slice', 'Fit X-slice', 'Fit Y-slice'}, 'Location', 'northeast' , 'FontSize', 10)
     set(gcf,'Name',['frame = ',num2str(frame), ' spot # ',num2str(ii)]) 
     title ('Gauss + erf fit', 'FontSize', 10)

    disp('... paused... Click any button to proceded to the next spot ') 
    pause
end

    
    new_spotData=spotData;
    
    params.script=mfilename;
 params.date=date;
 params.version=version;
params.fitRadius = fitRadius;
params.marginRadius = marginRadius;
params.imageScale = imageScale;
params.markerSize = markerSize;
params.fitWeight = fitWeight;
params.showFitData = showFitData;
    
end
% 
% spotFinder spotList data'
%             Spots=spotList{frame};
%             spotData=Spots{spotID};
%             h=Spots{spotID}.h;
%             w=Spots{spotID}.w;
%             w2=w^2;
%             b=Spots{spotID}.b;
%             x=Spots{spotID}.x;
%             y=Spots{spotID}.y;
%             pos=round([x,y]);




