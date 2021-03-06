%% riv_pantanalMAP Demo - This demo is a walkthrough demonstrating many of the
% functions included in the riv_pantanalMAP toolbox. The demo is best experienced by
% running blocks of code at a time by highlighting the lines to run and
% either pressing F9 or right-clicking and choosing "Evaluate Selection."
% A supplemental .docx file is also provided to show results and provide
% further explanation. Cell numbering in this script corresonds to
% numbering in that document.

%% 1 - Load data
% Load the data to analyze - The data are from the middle of R6. See the 
% reach in Google Earth Engine: https://earthengine.google.com/timelapse/#v=-9.5,-74.13468,9.148,latLng&t=2.73
% Specifically, a bounding box was drawn in R6 with the following limits:
% xl = [706 2050];
% yl = [1715 3537];
% The entire R6 box is not provided in this demo. You must change the path
% to the riv_pantanalMAP folder you downloaded.
RivMAP_path = 'C:\Users\renan\Unifesp\Iniciação Científica\ESTIMATIVA DA CONCENTRAÇÃO DE SEDIMENTOS EM SUSPENSÃO NA BACIA HIDROGRÁFICA DO RIO TAQUARIMS\RivMAP';
cd(RivMAP_path)

load('pantanal/riv_pantanal.mat')

% Clear the workspace
close all; clearvars -except riv_pantanal; clc;

% Turn off image size plotting warning
warning('off','images:initSize:adjustingMag');

% Data are stored in a structure called 'riv_pantanal.' There are currently two
% substructures in riv_pantanal: 'meta' and 'im.' meta contains the year, exit
% sides, and nominal channel width of the channel. Run the following lines
% to see what these values are for 1984.
disp(['Exit sides are ',num2str(riv_pantanal(1).meta.exit_sides),'.']);
disp(['Year is ',num2str(riv_pantanal(1).meta.year),'.']);
disp(['Nominal channel width is ',num2str(riv_pantanal(1).meta.Wn),'.']);

% In the 'im' substructure, there are two binary channel masks: 'hc'
% (hydraulically-connected) and 'st' (single-thread). The 'st' images were
% created from the 'hc' images by removing secondary channels.

%% 2. Plot channel masks

% Let's look at the 32 years of hydraulically-connected channel mask stored
% in the 'im.hc' substructure.
close all
for i = 1:numel(riv_pantanal)
    I = riv_pantanal(i).im.st;
    imshow(I)
    title(riv_pantanal(i).meta.year)
    pause(0.1)
end

%% 4. Process all years
% If you want to skip this processing, you can load in the variable
% 'riv_pantanal_processed' contained in the riv_pantanalMAP demo folder.
plotornot = 0; % Don't plot
for i = 1:numel(riv_pantanal)
    i
    % Load variables
    Wn = riv_pantanal(i).meta.Wn;
    I = riv_pantanal(i).im.st;
    es = riv_pantanal(i).meta.exit_sides;

    % Centerline
    [cl, Icl] = centerline_from_mask(I,es,Wn,plotornot);
    cls = savfilt(cl,Wn);

    % Banklines
    banks = banklines_from_mask(I, riv_pantanal(i).meta.exit_sides, plotornot);
    lb = banks{1}; 
    rb = banks{2};
    
    % Smooth banklines
    lbs = savfilt(lb,round(1.5*Wn));
    rbs = savfilt(rb,round(1.5*Wn));

    % Channel width
    W = width_from_banklines(cl, lb, rb, Wn);
    Wavg = mean(W(isnan(W)==0)); % average of nodewise widths

    % Distance between nodes (along-stream distance)
    S = [0; cumsum(sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2))];

    % Channel directions
    A = angles(cls);

    % Centerline curvatures
    C = curvatures(cls);

    % Reach average width
    Wra = sum(sum(I))/S(end);

    % Save variables - we create a new substructure called 'vec' to store
    % vector-based variables
    riv_pantanal(i).vec.cl = cl;
    riv_pantanal(i).vec.cls = cls;
    riv_pantanal(i).vec.lb = lb;
    riv_pantanal(i).vec.lbs = lbs;
    riv_pantanal(i).vec.rb = rb;
    riv_pantanal(i).vec.rbs = rbs;

    riv_pantanal(i).vec.W = W;
    riv_pantanal(i).vec.Wavg = Wavg; % Average of pointwise widths
    riv_pantanal(i).vec.Wra = Wra; % Reach average width
    riv_pantanal(i).vec.S = S;
    riv_pantanal(i).vec.A = A;
    riv_pantanal(i).vec.C = C;
    riv_pantanal(i).vec.cl_len = S(end); % Total centerline length

    riv_pantanal(i).im.cl = Icl; % Store the image of the centerline also

    disp(['Year ',num2str(i),'/',num2str(numel(riv_pantanal)),' is finished.'])
end
% Clear the workspace
close all; clearvars -except riv_pantanal; clc;
save('riv_pantanal','riv_pantanal')

%% 3. Analyze Centerline From Mask

% Now we will walk through one year of planform analysis
plotornot = 1; % We want to plot our results
close all
for i = 1:numel(riv_pantanal)
    % Load the variables we need from the riv_pantanal structure
    Wn = riv_pantanal(i).meta.Wn; % nominal width
    Ist = riv_pantanal(i).im.st; % single-thread channel mask
    es = riv_pantanal(i).meta.exit_sides; % exit sides

    centerline_from_mask(Ist,es,Wn,plotornot)
    title(riv_pantanal(i).meta.year)
    pause(0.1)
end

% Clear the workspace
close all; clearvars -except riv_pantanal; clc;

%% 4. Analyze Banklines From Mask

% Now we will walk through one year of planform analysis
plotornot = 1; % We want to plot our results
close all
for i = 1:numel(riv_pantanal)
    
    % Load the variables we need from the riv_pantanal structure
    Ist = riv_pantanal(i).im.st; % single-thread channel mask
    es = riv_pantanal(i).meta.exit_sides; % exit sides

    banklines_from_mask(Ist,es,plotornot)
    title(riv_pantanal(i).meta.year)
    pause(0.1)
end

% Clear the workspace
close all; clearvars -except riv_pantanal; clc;

%% 5. Analyze Angles

% Now we will walk through one year of planform analysis
plotornot = 1; % We want to plot our results
close all
for i = 1:numel(riv_pantanal)
    disp(i)
     Wn = riv_pantanal(i).meta.Wn; % nominal width
     cl = riv_pantanal(i).vec.cl;
     cls = savfilt(cl,Wn);

     % Let's plot the two width methods to see how they compare
     % First, we need to compute the streamwise distance along the centerline
     S = [0; cumsum(sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2))];

     % Now plot both widths
     close all
     %plot(S,width_banklines(i).wbl); hold on
     %plot(width_from_mask(i).Icl,width_from_mask(i).cl,'r');
     xlabel('streamwise distance, pixels'); ylabel('width, pixels')
     legend('W_c_l','W_m_a_s_k')

     % Channel directions
     A = angles(cls);

     % Centerline curvatures
     C = curvatures(cls);

     % Plot angles and curvatures
     close all;
     subplot(2,1,1)
     plot(S,A);
     title(riv_pantanal(i).meta.year)
     xlabel('streamwise distance, pixels'); ylabel('channel direction, radians')
     subplot(2,1,2)
     plot(S,C,'r');
     ylim([-0.2 0.2])
     xlabel('streamwise distance, pixels'); ylabel('curvature, pixels^-^1')
     pause(2)
end

% Clear the workspace
close all; clearvars -except riv_pantanal; clc;

%% 6. Analyze Migration
for i = 1:numel(riv_pantanal)-1 %
    
    % Find index of next time step (for cases where data is missing from a
    % year)
    for jj = (i+1):numel(riv_pantanal)
        if isempty(riv_pantanal(jj).vec) == 0
            i2 = jj;
            break
        end
    end
    
    % Compute centerline migration areas and find cutoffs
    cl1 = riv_pantanal(i).vec.cls; % We are using the smooth centerlines but can also use unsmoothed
    cl2 = riv_pantanal(i2).vec.cls;
    Wn = riv_pantanal(i).meta.Wn;
    es1 = riv_pantanal(i).meta.exit_sides;
    es2 = riv_pantanal(i2).meta.exit_sides;
    sizeI = size(riv_pantanal(i).im.cl);
    [Imig,Icutsall,cutidcs,cutarea,cutlen,chutelen] = migration_cl(cl1, cl2, es1, es2, Wn, sizeI);
    % Compute erosion/accretion areas and associated cutoffs
    I1 = riv_pantanal(i).im.st;
    I2 = riv_pantanal(i2).im.st;
    [Ie,Ia,Inc,Icutsm] = migration_mask(I1,I2,Wn);
    
    % Store results
    riv_pantanal(i).mig.cl.Imig = Imig;
    riv_pantanal(i).mig.cl.Icuts = Icutsall;
    riv_pantanal(i).mig.cl.cutidcs = cutidcs;
    riv_pantanal(i).mig.cl.cutareas = cutarea;
    riv_pantanal(i).mig.cl.cutlen = cutlen;
    riv_pantanal(i).mig.cl.chutelen = chutelen;

    
    riv_pantanal(i).mig.mask.Ie = Ie;
    riv_pantanal(i).mig.mask.Ia = Ia;
    riv_pantanal(i).mig.mask.Inc = Inc;
    riv_pantanal(i).mig.mask.Icuts = Icutsm;
    
    disp(['Year ',num2str(i),'/',num2str(numel(riv_pantanal)-1),' is finished.'])
end

%% 7. Plot migration maps

% At this point, we have computed all the planform changes supported by 
% riv_pantanalMAP. Now we will look at some results.

% First, let's look at the total migrated area across all times. We will
% initialize a blank image, loop through all the migrated images and mark
% each migrated pixel as true (1) in the blank image.
Imig = false(size(riv_pantanal(1).im.st)); % initialize image
Ie = false(size(riv_pantanal(1).im.st)); % initialize image
Ia = false(size(riv_pantanal(1).im.st)); % initialize image
Icuts = false(size(riv_pantanal(1).im.st)); % initialize image
for i = 1:numel(riv_pantanal)-1
    Imig(riv_pantanal(i).mig.cl.Imig) = true;
    Ie(riv_pantanal(i).mig.mask.Ie) = true;
    Ia(riv_pantanal(i).mig.mask.Ia) = true;
    Icuts(riv_pantanal(i).mig.cl.Icuts) = true;
    yrs(i) = riv_pantanal(i).meta.year;
end

% Centerline migrated area map
close all
subplot(1,3,1)
subimage(Imig)
title('Migrated area (cl)','fontsize',16)

% We can also plot the areas colored as a function of time
Imigc = zeros(size(riv_pantanal(1).im.st)); % initialize image (use zeros instead of false so that we can store non-binary numbers)
for i = 1:numel(riv_pantanal)-1
    Imigc(riv_pantanal(i).mig.cl.Imig) = i;
    yrs(i) = riv_pantanal(i).meta.year;
end
subplot(1,3,2)
cmap = colormap(parula(numel(riv_pantanal)));
subimage(Imigc, cmap);
cb = colorbar;
ytix = [2:5:30]/(numel(riv_pantanal)-1);
ytixl = yrs(1:3:6);
set(cb,'ytick',ytix,'yticklabel',ytixl)
title('Migrated area, colored by year','fontsize',16)

% Erosion/accretion map
subplot(1,3,3)
Iaeplot = imfuse(Ie,Ia);
subimage(Iaeplot)
title('Erosion (green) and accretion (magenta)','fontsize',16)
set(gcf, 'Position', get(0,'Screensize')-[-45 -30 120 120]); % Maximize figure size

%% 8. Plot cutoffs
close all
subplot(1,3,1);
subimage(Icuts) % It appears that five cutoffs have occurred
title('Binary map of cutoffs','fontsize',16)

% Let's plot the cutoffs again, but color them according to year
Icutc = zeros(size(riv_pantanal(1).im.st)); % initialize image (use zeros instead of false so that we can store non-binary numbers)
for i = 1:numel(riv_pantanal)-1
    Icutc(riv_pantanal(i).mig.cl.Icuts) = i;
    yrs(i) = riv_pantanal(i).meta.year;
end
subplot(1,3,2);
subimage(Icutc,cmap)
cb = colorbar;
ytix = [2:5:30]/(numel(riv_pantanal)-1);
ytixl = yrs(1:3:6);
set(cb,'ytick',ytix,'yticklabel',ytixl)
title('Migrated area, colored by year','fontsize',16)
title('Cutoffs colored by year','fontsize',16)

% Now it appears there were six cutoffs. Let's plot cutoff areas in time to
% make sure.
for i = 1:numel(riv_pantanal)-1
    n = numel(riv_pantanal(i).mig.cl.cutareas);
    if n > 0
        cutareas(i,:) = riv_pantanal(i).mig.cl.cutareas;
    end    
end
subplot(1,3,3)
b1 = bar(yrs,cutareas*900/10^6);
b1.FaceColor = 'b';
xlim([1983 2015])
ylabel('Cutoff area, km^2')
title('Cutoff areas through time')
set(gca,'fontsize',16)
set(gcf, 'Position', get(0,'Screensize')-[-45 -30 120 120]); % Maximize figure size
% Now we see there were seven total cutoffs.
