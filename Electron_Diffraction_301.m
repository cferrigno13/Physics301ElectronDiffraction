%% Physics 301 Electron Diffraction Lab


%% Importing Raw Data

Vt = [1.5:.25:5.5];
V = 1000 * Vt';

Ds = [36.71;34.85;30.33;29.07;28.33;26.27;25.65;24.83;22.48;22.73;21.74;21.13;20.77;20.29;20.05;19.49;18.39];
Db = [63.82;58.01;54.48;51.64;49.69;46.93;45.38;43.05;42.77;40.81;39.10;38.26;36.66;35.64;34.79;33.54;33.07];

L = 130; %+/- 2mm
R = 66; %+/- .1mm
tt = 1.5;

%% Plotting for the small ring

thetaS = 0.5 * atan(Ds./(2*(L-R+(R^2-(Ds/2).^2).^(1/2))));

%% Find Uncertainty in Theta

%Uncertainty due to R uncertainty
Rs = R + .1;
UfromR = thetaS - 0.5 * atan(Ds./(2*(L-Rs+(Rs^2-(Ds/2).^2).^(1/2))));

%Uncertainty due to L uncertainty
Ls = L + 2;
UfromL = thetaS - 0.5 * atan(Ds./(2*(Ls-R+(R^2-(Ds/2).^2).^(1/2))));

%Uncertainty due to D uncertainty
Dss = Ds + 2;
UfromDs = thetaS - 0.5 * atan(Dss./(2*(L-R+(R^2-(Dss/2).^2).^(1/2))));

%Add in quadrature
thetaSerror = (UfromR.^2 + UfromL.^2 + UfromDs.^2).^(1/2);


%% Variables to Plot 

Ys = sin(thetaS);
Xs = V.^(-1/2);

%% Uncertainty in Ys

UYs = Ys - sin(thetaS + thetaSerror);

%% Plotting and Fitting

ps = polyfit(Xs,Ys,1);
PS = Xs.*ps(1) + ps(2);
d2 = .614 / ps(1);

% Find Chi Squared
chi2s = 0;
for i = 1:length(Xs)
    chi2s = chi2s + ((Ys(i) - PS(i))^2 / (UYs(i)^2));
end

LittleRing = figure(1);

errorbar(Xs,Ys,UYs,'b-o')
hold on
plot(Xs,ps(1)*Xs+ps(2))
ylabel('sin(theta)')
xlabel('Voltage^-^1^/^2')
title('Small Ring: Plot to find d2')
legend('Observed Data','Linear Fit')
hold off
xl = xlim;
yl = ylim;
xt = .05 * (xl(2)-xl(1)) + xl(1);
yt = .9 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x + %f \nd2 = {%f} nm \nChi^2 = {%f}',ps(1),ps(2),d2,chi2s);
text(xt,yt,caption,'FontSize',14,'Color', 'k', 'FontWeight', 'bold');




saveas(LittleRing,'SmallRing_Plotd2.png')

%% Plotting for the big ring

thetaB = 0.5 * atan(Db./(2*(L-R+(R^2-(Db/2).^2).^(1/2))));
%% Find Uncertainty in Theta

%Uncertainty due to R uncertainty
Rs = R + .1;
UfromR = thetaB - 0.5 * atan(Db./(2*(L-Rs+(Rs^2-(Db/2).^2).^(1/2))));

%Uncertainty due to L uncertainty
Ls = L + 2;
UfromL = thetaB - 0.5 * atan(Db./(2*(Ls-R+(R^2-(Db/2).^2).^(1/2))));

%Uncertainty due to D uncertainty
Dbs = Db + 2;
UfromDb = thetaB - 0.5 * atan(Dbs./(2*(L-R+(R^2-(Dbs/2).^2).^(1/2))));

%Add in quadrature
thetaBerror = (UfromR.^2 + UfromL.^2 + UfromDb.^2).^(1/2);

%% Variables to Plot
Yb = sin(thetaB);
Xb = V.^(-1/2);

%% Uncertainty in Yb

UYb = Yb - sin(thetaB + thetaBerror);

%% Plotting and Fitting

pb = polyfit(Xb,Yb,1);
PB = Xb.*pb(1) + pb(2);
d1 = .614 / pb(1);

% Find Chi Squared
chi2b = 0;
for i = 1:length(Xb)
    chi2b = chi2b + ((Yb(i) - PB(i))^2 / (UYb(i)^2));
end

BigRing = figure(2);

errorbar(Xb,Yb,UYb,'b-o')
hold on
plot(Xb,pb(1)*Xb+pb(2))
ylabel('sin(theta)')
xlabel('Voltage^-^1^/^2')
title('Big Ring: Plot to find d1')
legend('Observed Data','Linear Fit')
hold off
xl = xlim;
yl = ylim;
xt = .05 * (xl(2)-xl(1)) + xl(1);
yt = .9 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x + %f \nd1 = {%f} nm \nChi^2 = {%f}',pb(1),pb(2),d1,chi2b);
text(xt,yt,caption,'FontSize',14,'Color', 'k', 'FontWeight', 'bold');

saveas(BigRing,'BigRing_Plotd1.png')
