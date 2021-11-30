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

%% 3. Analyze one year of planform characteristics

% Now we will walk through one year of planform analysis
i = 1; % We will analyze the first element in the structure (can change i if you want to analyze a different year)
plotornot = 1; % We want to plot our results

% Load the variables we need from the riv_pantanal structure
Wn = riv_pantanal(i).meta.Wn; % nominal width
Ist = riv_pantanal(i).im.st; % single-thread channel mask
Ihc = riv_pantanal(i).im.hc; % hydraulically-connected channel mask
es = riv_pantanal(i).meta.exit_sides; % exit sides

%% 3a. Compute the centerline
close all
[cl, Icl] = centerline_from_mask(Ist,es,Wn,plotornot);

%% 3b. Compute the banklines
close all
banks = banklines_from_mask(Ist, es, plotornot);

%% 3c. Along-channel widths

% Channel width at each pixel from banklines
Wbl = width_from_banklines(cl, banks{1}, banks{2}, Wn);

Wavg = mean(Wbl(isnan(Wbl)==0)); % average of pixel-wise widths - note that it doesn't quite agree with our nominal channel width (Wn), but this is fine since Wn is just used to parameterize buffer boxes

% Channel width from mask at specified spacing
spacing = Wn/2;
[Wm, SWm] = width_from_mask(Ist, cl, spacing);

% Let's plot the two width methods to see how they compare
% First, we need to compute the streamwise distance along the centerline
S = [0; cumsum(sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2))];

% Now plot both widths
close all
plot(S,Wbl); hold on
plot(SWm,Wm,'r');
xlabel('streamwise distance, pixels'); ylabel('width, pixels')
legend('W_c_l','W_m_a_s_k')

%% 3d. Angles and curvatures

% Now we will compute channel direction (angles) and curvatures of the
% centerline. We will first smooth the centerline.
cls = savfilt(cl,Wn);

% Channel directions
A = angles(cls);

% Centerline curvatures
C = curvatures(cls);

% Plot angles and curvatures
close all;
subplot(2,1,1)
plot(S,A);
xlabel('streamwise distance, pixels'); ylabel('channel direction, radians')
subplot(2,1,2)
plot(S,C,'r');
ylim([-0.2 0.2])
xlabel('streamwise distance, pixels'); ylabel('curvature, pixels^-^1')

% Two things of note: (1) Curvature is quite noisy. This is typical when
% estimating second-deriv_pantanalative quantities. For smoother curvatures, the
% centerline could be smoothed mulitple times and/or with a larger window.
% (2) The ends of both signals are somewhat noisy--this is an artifact of
% how the savfilt smoother handles end cases. 

%% 3e. Reach-averaged widths

% Finally, we can compute the average width for the entire reach
Wra = sum(sum(Ist))/S(end);

% This value agrees fairly well with the average of widths computed on a
% pixel-wise basis:
Wra
Wavg

%% 4. Process all years
% If you want to skip this processing, you can load in the variable
% 'riv_pantanal_processed' contained in the riv_pantanalMAP demo folder.
plotornot = 0; % Don't plot
for i = 1:numel(riv_pantanal)
    
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
save('riv_pantanal','riv_pantanal')

%% 5. Plot all years
close all
for i = 1:numel(riv_pantanal)
    imshow(riv_pantanal(i).im.st); hold on
    plot(riv_pantanal(i).vec.cls(:,1),riv_pantanal(i).vec.cls(:,2),'b','linewidth',1.5);
    plot(riv_pantanal(i).vec.lb(:,1),riv_pantanal(i).vec.lb(:,2),'m','linewidth',1.5)
    plot(riv_pantanal(i).vec.rb(:,1),riv_pantanal(i).vec.rb(:,2),'m','linewidth',1.5)
    title(riv_pantanal(i).meta.year)
    pause(0.2)
end

% Note that if you want to plot banklines, centerlines without plotting on
% top of an image, you need to multiply the y-coordinates by -1. This is
% because the origin for images is at the top-left of the image, while for
% ordinary plotting the origin is at the bottom-left. E.g.

% close all
% plot(riv_pantanal(1).vec.cls(:,1),-riv_pantanal(1).vec.cls(:,2),'b','linewidth',1.5); axis equal; hold on
% plot(riv_pantanal(1).vec.lb(:,1),-riv_pantanal(1).vec.lb(:,2),'m','linewidth',1.5)
% plot(riv_pantanal(1).vec.rb(:,1),-riv_pantanal(1).vec.rb(:,2),'m','linewidth',1.5)


%% 6. Compute migrations

% Compute changes in time (migrations, erosions, accrections, cutoffs).
% Results related to migration are stored in a new substructure called
% 'mig.' This substructure contains two sub-substructures: (1) 'cl' is used
% to store results related to centerline migrations, and (2) 'mask' stores
% results from mask differencing (erosion/accretion).
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
save('riv_pantanal','riv_pantanal')

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
ytixl = yrs(2:5:30);
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
ytixl = yrs(2:5:30);
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


%% 9. Annual reachwide migration rates

% Let's start by computing the average migration rate for the entire reach
% across all time. Average migration rate is computed as the migrated area
% divided by the centerline length--we have both stored in our riv_pantanal 
% structure.
for i = 1:numel(riv_pantanal)-1
    MAcl(i) = sum(sum(riv_pantanal(i).mig.cl.Imig));
    MAe(i) = sum(sum(riv_pantanal(i).mig.mask.Ie));
    MAa(i) = sum(sum(riv_pantanal(i).mig.mask.Ia));
    len(i) = riv_pantanal(i).vec.cl_len;
end
Mrcl = MAcl./len;
Mre = MAe./len;
Mra = MAa./len;

%% 9a. Reachwide migration through time
close all
plot(yrs,Mrcl*30,'r'); hold on
plot(yrs,Mre*30,'b');
legend('CL','Erosion','Accretion','location','best')
ylabel('Reach average migration rate, m/yr')
xlim([1984 2015])
set(gca,'fontsize',16)

%% 9b. Spatial variation of migration rates

% We can also see how migration rates change spatially by defining a
% meander-belt centerline and averaging migrations through time along this
% centerline. 

% The riv_pantanalMAP function spatial_migration will compute a meander belt that
% encompasses all the migrated area, as well as a centerline of the meander
% belt. It returns the average migration rate along the meander belt
% centerline. We need to feed it a binary image of all channel positions
% along with the image of all migrations over the time period, and a cell
% containing all the centerline images.
clear Icp Icl Imig Ie Ia
for i = 1:numel(riv_pantanal)
    Icp{i} = riv_pantanal(i).im.st;
    Icl{i} = riv_pantanal(i).im.cl;
end
for i = 1:numel(riv_pantanal)-1
   Imig{i} = riv_pantanal(i).mig.cl.Imig;
   Ie{i} = riv_pantanal(i).mig.mask.Ie;
   Ia{i} = riv_pantanal(i).mig.mask.Ia;
end

Wn = riv_pantanal(1).meta.Wn;
es = riv_pantanal(1).meta.exit_sides;
spacing = 2.1*Wn; % we'll compute migration rates every two channel widths
plotornot = 0;

% Populate our input cell for analyzing with spatial changes
Ianalyze{1} = Imig;
Ianalyze{2} = Ie;
Ianalyze{3} = Ia;

% Now we will call spatial_changes.
[Iout, Achan, cllen, belt] = spatial_changes(Icp, Icl, Ianalyze, spacing, Wn, es, plotornot);
% Unpack the results
Acl = Iout{1};
AE = Iout{2};
AA = Iout{3};
% To compute annual average migration rates, we divide Atot by lenavg, then 
% divide by the total number of elapsed years. Multiply by 30 to convert
% from pixels to meters.
nyrs = numel(riv_pantanal)-1; % number of years of migrated area in the Imig image

% The variable 'belt' contains the info we need to reconstruct the meander
% belt boundaries and buffer polygons. Let's see what the segments look 
% like.
% First, make an image of all channel positions through time
Icp_all = false(size(Icp{1}));
Imig_all = Icp_all;
for i = 1:numel(Imig)
    Icp_all(Icp{i}) = true;
    Imig_all(Imig{i}) = true;
end
% Now plot the buffer polygons
close all
subplot(1,3,1)
subimage(Icp_all); hold on
plot(belt.lb(:,1),belt.lb(:,2),'m');
plot(belt.rb(:,1),belt.rb(:,2),'m');
for i = 1:length(belt.cl)
    plot([belt.lb(i,1) belt.cl(i,1) belt.rb(i,1)],[belt.lb(i,2) belt.cl(i,2) belt.rb(i,2)],'m')
end
xlim([0 size(Icp_all,2)])
ylim([0 size(Icp_all,1)])
title('Meander belt and segments','fontsize',16)

% Let's plot the results
subplot(1,3,2)
subimage(Imig_all); hold on
plot(belt.cl(:,1),belt.cl(:,2),'r.','markersize',8)
for i = 1:length(belt.cl)
    text(belt.cl(i,1)+5,belt.cl(i,2)+5,num2str(i),'color','r')
end
title('Meander belt centerline nodes','fontsize',16)

% Now plot the along-channel migration rates
M = sum(Acl./cllen,1)/nyrs*30;

subplot(1,3,3)
plot(M,1:numel(M),'k')
ylim([0 numel(M)])
ylabel('Along-stream distance in widths')
xlabel('Average annual migration rate, m/yr')
set(gca,'fontsize',14)
title('Spatial variation in migration rate','fontsize',16)
set(gcf, 'Position', get(0,'Screensize')-[-45 -30 120 120]); % Maximize figure size

%% 9c. Erosion and accretion rates

% We can also look at how erosion and accretion change in tandem
close all
subplot(3,1,1)
MAe = sum(AE./cllen,1)*30; % Conversion to m from pixels
MAa = sum(AA./cllen,1)*30;
plot(yrs,MAe,'k'); hold on 
plot(yrs,MAa,'r');
ylabel('Rate, m/yr')
axis tight
legend('Erosion','Accretion','location','best');
set(gca,'fontsize',14)
title('Erosion and accretion rates','fontsize',16)
xlim([1984 2014])

% We can see the erosion/accretion balance better by looking at their
% cumulative values for the reach
subplot(3,1,2)
plot(yrs,cumsum(MAe)*900/10^6,'k'); hold on % Converted to km^2
plot(yrs,cumsum(MAa)*900/10^6,'r')
legend('Erosion','Accretion','location','best')
ylabel('Cumulative area, km^2')
set(gca,'fontsize',14)
title('Cumulative erosion and accretion','fontsize',16)
xlim([1984 2014])

% Over this time period, nearly 100 km^2 more of erosion occurred than
% accretion. This indicates channel widening from 1996 onward, so we should
% see this if we plot average channel width. We will also plot the
% difference between the cumulative sums of erosion and accrection.
for i = 1:numel(riv_pantanal)
    Wr_avg(i) = riv_pantanal(i).vec.Wavg;
end
subplot(3,1,3)
[hAx,hLine1,hLine2] = plotyy(yrs,Wr_avg(1:31),yrs,(cumsum(MAe)-cumsum(MAa))*900/10^6);
hLine2.LineStyle = ':';
hLine2.LineWidth = 2;
hLine1.Color = 'k';
hLine2.Color = 'b';
hAx(1).YColor = 'k';
hAx(2).YColor = 'b';
set(hAx,'fontsize',14);
ylabel(hAx(2),'Cumulative(E) - Cumulative(A), km^2') % right y-axis
ylabel(hAx(1),'Average width, pixels') % left y-axis
ylim(hAx(2),[-10 25])
ylim(hAx(1),[18 30])
xlim([1984 2014])
leg = legend('Width','E-A','location','best');
set(leg,'fontsize',16)
title('Width and net erosion/accretion','fontsize',16)
set(gcf, 'Position', [46    31   729   960]); 

%% 10. Spacetime maps
% Make an image of areas computed by spatial_changes. The x-axis will be
% year; the y-axis will be the jth buffer polygon. We can compute the
% along-centerline distance of each buffer polygon from the belt structure
% and re-label the y-axis with along-stream distance.
close all
subplot(1,2,1)
map1 = bone(1500);
h = subimage(Achan./cllen*30,map1);
caxis([500 1500])
S = belt.Smid*30/1000;
xtixl = 1985:5:2014;
xtix = xtixl-min(xtixl)+1;
caxis([300 1000])
set(gca,'YDir','normal','xtick',xtix,'xticklabel',xtixl,'yticklabel',S)
ylabel('distance along belt centerline , km')
title('width, m')
set(gca,'fontsize',14)

subplot(1,2,2)
map2 = parula(30);
h = subimage(Acl./Achan*100,map2);
set(gca,'YDir','normal','xtick',xtix,'xticklabel',xtixl,'yticklabel',S)
ylabel('distance along belt centerline , km')
title('centerline migration rate normalized by channel area, %/yr')

set(gcf, 'Position', get(0,'Screensize')-[-45 -30 120 120]); % Maximize figure size
set(gca,'fontsize',14)

%% 11. Planform changes for a portion of the reach
% If we are interested in only a smaller sub-reach within our larger reach,
% we can define a bounding box around the sub-reach and perform the
% analysis only within this box.

% We will focus on the downstream portion of the reach.
xl = [150 600];
yl = [14 560];
I = riv_pantanal(1).im.st;
Imask = false(size(I));
Imask(yl(1):yl(2),xl(1):xl(2)) = true;
Icrop = I & Imask;
% close all
% imshowpair(I,Icrop)
% title('White indicates sub-reach of interest')

% We can use the migrated images we've already stored. We just need to
% recompute the centerline length for the subreach so we can calculate
% migration rates. We will use the riv_pantanal.im.cl image we already computed to
% find the subreach centerline length.
clear MAcl MAe MAa
for i = 1:numel(riv_pantanal)-1
    % First find the subreach centerline length
    Icl = riv_pantanal(i).im.cl; % load the centerline image
    Icl(~Imask) = false; % mask the centerline image
    E = find(bwmorph(Icl,'endpoints')); % find endpoints of the centerline
    D = bwdistgeodesic(Icl,E(1),'quasi'); % compute euclidean distances along centerline from endpoint
    subreach_len(i) = max(max(D)); % centerline length is max of distances (other endpoint)
    
    % Now find migration areas within the subreach
    Imig2 = riv_pantanal(i).mig.cl.Imig;
    Ie = riv_pantanal(i).mig.mask.Ie;
    Ia = riv_pantanal(i).mig.mask.Ia;
    MAcl(i) = sum(sum(Imig2(Imask)));
    MAe(i) = sum(sum(Ie(Imask)));
    MAa(i) = sum(sum(Ia(Imask))); 
end
Mrcl_sr = MAcl./subreach_len; % centerline migration rate
Mre_sr = MAe./subreach_len; % erosion rate
Mra_sr = MAa./subreach_len; % accretion rate

close all
subplot(3,1,1)
plot(yrs,Mrcl,'r'); hold on
plot(yrs,Mrcl_sr,'k');
xlim([1984 2015])
ylabel('Migration rate (CL), m/yr','fontsize',16)
leg = legend('Sub-reach','Entire reach','location','best');
set(leg,'fontsize',16)
set(gca,'fontsize',14)
title('Migration rate for sub-reach and entire reach')

% From 2005-2010, the sub-reach experienced significantly larger migration
% rates than the entire reach. Could this be due to a cutoff upstream of
% the reach? From our cutoffs-in-time plot, we saw there was a cutoff in
% 2004. Let's see if it was close to our sub-reach.
idx = find(yrs==2004);
Icut2004 = (riv_pantanal(idx).mig.cl.Icuts);
I = riv_pantanal(idx).im.st;
Isubreach = I & Imask;
I2004 = zeros(size(I));
I2004(I) = 1;
I2004(Isubreach) = 2;
I2004(Icut2004) = 3;
subplot(3,1,[2 3])
imshow(I2004)
cmap = [0 0 0; 1 1 1; 1 0 0; 0 0 1];
colormap(cmap)
caxis([0 3])
title('2004: Red = subreach, Blue = cutoff','fontsize',16)
set(gcf, 'Position', get(0,'Screensize')-[-45 -30 120 120]); % Maximize figure size

% The cutoff in 2004 was directly upstream of our subreach. It is possible
% that this cutoff induced the accelerated migration rates within the
% subreach that we observed from 2005-2010.

%% 12. Combining georeferenced images
% Load in the a separate dataset that contains imagery from four large
% boxes of the Ucayali riv_pantanaler
load('riv_pantanal_georeffed')
clear I G
% We need to create two cell variables: one contains each of the images we
% want to combine, and one contains the images' georeferencing info.
for i = 1:4
    if i == 1
        I{i} = riv_pantanal3.im.st;
        G{i} = riv_pantanal3.meta.georef;
    elseif i == 2
        I{i} = riv_pantanal4.im.st;
        G{i} = riv_pantanal4.meta.georef;
    elseif i == 3
        I{i} = riv_pantanal5.im.st;
        G{i} = riv_pantanal5.meta.georef;
    elseif i == 4
        I{i} = riv_pantanal6.im.st;
        G{i} = riv_pantanal6.meta.georef;
    end
end
% Stitch them together
[Im_out,Georef_out] = combine_georeffed_images(I,G);

% Plot the individual images and their stitched versions
close all
figure;
title('Channel mask images')
subplot(2,3,1)
subimage(I{1})
title('R3')
subplot(2,3,2)
subimage(I{2})
title('R4')
subplot(2,3,4)
subimage(I{3})
title('R5')
subplot(2,3,5)
subimage(I{4})
title('R6')
subplot(2,3,[3 6])
subimage(Im_out)
title('All')

% We can also combine centerlines in vector format using
% combine_georeffed_images
clear I G
for i = 1:4
    if i == 1
        I{i} = riv_pantanal3.vec.cls;
        G{i} = riv_pantanal3.meta.georef;
    elseif i == 2
        I{i} = riv_pantanal4.vec.cls;
        G{i} = riv_pantanal4.meta.georef;
    elseif i == 3
        I{i} = riv_pantanal5.vec.cls;
        G{i} = riv_pantanal5.meta.georef;
    elseif i == 4
        I{i} = riv_pantanal6.vec.cls;
        G{i} = riv_pantanal6.meta.georef;
    end
end
[cl_out,Georef_out] = combine_georeffed_images(I,G);
figure
title('Centerlines')
subplot(2,3,1)
plot(I{1}(:,1),-I{1}(:,2),'b'); axis equal;
title('R3')
subplot(2,3,2)
plot(I{2}(:,1),-I{2}(:,2),'b'); axis equal;
title('R4')
subplot(2,3,4)
plot(I{3}(:,1),-I{3}(:,2),'b'); axis equal;
title('R5')
subplot(2,3,5)
plot(I{4}(:,1),-I{4}(:,2),'b'); axis equal;
title('R6')
subplot(2,3,[3 6])
plot(cl_out(:,1),-cl_out(:,2),'b'); axis equal;
title('All')

