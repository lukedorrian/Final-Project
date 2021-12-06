load('Atlanta_O3.txt')
Atlanta_O3 = Atlanta_O3(Atlanta_O3(:,3)>0,:); %filter out negative values
measurements = Atlanta_O3(:,3)';
simulated = Atlanta_O3(:,4)';

[p,s] = polyfit(measurements, simulated, 1);
slope = p(1);
intercept = p(2);
tValue = abs(tinv(0.025,length(measurements) - 2));
[y_fit, delta] = polyval(p, measurements, s);

errorVariance = sum((simulated - y_fit).^2)/(length(simulated)-2);
varX = var(measurements);
sigmaSlope = sqrt(errorVariance/varX);
%95% confidence interval for line slope
slopeInterval = [slope-tValue*sigmaSlope slope+tValue*sigmaSlope];

[r,p, rlow, rhigh] = corrcoef(measurements, simulated);
b1 = r(1,2);
corrInterval = [rlow(1,2) rhigh(1,2)];
%significance test
significant = p(1,2) < 0.05;

%plotting
figure
plot(measurements,simulated,'bs', measurements, y_fit, 'r-')
hold on
plot(measurements, y_fit + 2.*delta, 'g-', measurements, y_fit - 2.*delta, 'g-');
title('Measurements vs Simulated O^3')

%assuming my derivation is correct

beta = sum(measurements.*simulated)/sum(measurements.^2);
newX = 0:0.5:max(measurements);
regressionLine = newX.*beta;
plot(newX, regressionLine, 'm')
hold off
