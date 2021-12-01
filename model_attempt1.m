%% initializations
rho = 917; % ice density (kg/m^3)
g = 9.81; % gravitational acceleration (m/s^2)
m = 2; % constant = 1/2(n + 1), so ~ 2
A = 200*1000*1.2; % basal roughness, need to consult Robel
beta = 2.5; % bed slope

% set dx and dt 
dx = 200; % meters
dt = 1; % years, may need to edit to make stepping more stable
t = 0:dt:1000; % years

% space (x) and thickness (h) 
h_i = 3000; % initial thickness (meters)
glacier_length = 50000; % glacier length (meters) ... somewhere around this
all_x = 0:dx:glacier_length; % position vector (meters)
all_h = zeros(length(t), length(all_x)); % thickness vector (meters)
% this is an array with the rows being time, columns being h(x,t)
all_h(1,1) = h_i; % set inland thickness 
all_h(1,:) = h_i - all_x.*(3*10^-3);

accum_i = 0.1; % initial accumulation rate (meters/year)
u_i = 20; % initial ice velocity (meters/year)

% need u * dt/dx < 1, close to 1 preferably
[h_test, a_test] = meshgrid(100:100:4000,-1:0.1:5);
u_test = (((rho*g/A)^m).*(h_test.^m).*(sind(a_test).^m));
u_stable = u_test >= 0.8 & u_test <= 1.2;
figure;
surf(h_test,a_test,u_test); hold on;
h = scatter3(h_test(u_stable),a_test(u_stable),u_test(u_stable),'filled', 'r');
colorbar;
title('Possible u Values (red where u >= 0.8 & <= 1.2)');
xlabel('h: thickness');
ylabel('a: alpha'); 
zlabel('u: ice velocity');

%% intialize outputs 
% go through each function h(x,t)
% calculate a(x) then u(x)
% right now timestepping by 1 year
dqdtdx(1:length(all_x)) = 0;
u_array = zeros(length(t),length(all_x));
a_array = u_array; % alpha -> surface slope

%% step through model 
for k = 1:(length(t) - 1) % time step 
    % thickness
    h_copy = all_h(k,:); % current thickness vector 
    h_diff = h_copy(1) - h_copy; 
    
    % sfc slope alpha
    a_apply(1) = beta; %how to incorporate a_apply?
    a_apply(2:length(all_x)) = beta - h_diff(2:end)./all_x(2:end); %how to incorporate proper a(x)
    a_array(k,:) = a_apply;

    % ice velocity 
    u_apply = (((rho*g/A)^m).*(h_copy.^m).*(sind(a_apply).^m));
    u_array(k,:) = u_apply;
    dqdtdx(1) = u_apply(1)*h_copy(1);

    % find dt*dq/dx
    for i = 2:length(all_x) % space step
        dqdtdx(i) = (dt/dx)*(h_copy(i)*u_apply(i) - h_copy(i-1)*u_apply(i-1));
    % multiplied by -1 since u_i controls mass flow OUT rather than in
    end
    
    all_h(k+1,:) = h_copy + dqdtdx + accum_i; 
end

figure;
plot(all_x, all_h(1,:),'r');
hold on;
plot(all_x, all_h(251,:), 'b'); 
plot(all_x, all_h(501,:), 'g');
plot(all_x, all_h(751,:), 'k');
plot(all_x, all_h(1001,:), 'm');
legend('h(x) at year 0', 'h(x) at year 250', 'h(x) at year 500', 'h(x) at year 750', 'h(x) at year 1000');
title(['H(x) over Time, Accumulation Rate ' num2str(accum_i) ' m/year']);
xlabel('x (m)');
ylabel('Glacier Depth (m)');
