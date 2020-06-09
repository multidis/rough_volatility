%% Clear workspace and load code:
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%% Define a rough Heston model (and inspect the object)
% Define a model:
model = rHestonClass('H',0.1,'nu',0.3,'rho',-0.65,'xi',0.15^2,'s0',100);

model

% Curves are stored as CurveClass objects:
model.y
model.q
model.xi

% Note how the yield curve (model.y) and dividend yield curve (model.q) by default 
% are assumed flat at zero. The forward variance curve (model.xi) is set to the fixed 
% value previously inputted. 

%% Compute implied volatility:
% Log-moneyness values, i.e. log(strike/forward):
k = (-0.3:0.05:0.2)'; 

% Expiries:
T = [0.05;0.2;0.5;1;2;3]; 

% Return prices for the cartesian product of the inputs:
cartProd = true;

% Compute prices:
ivSurf = model.GetPrices(k,T,cartProd);
ivSurf

% Plot:
for i=1:size(ivSurf.k,2)
    plot(ivSurf.k(:,i),ivSurf.iv(:,i),'-o','DisplayName',...
         ['T = ',num2str(T(i))],'LineWidth',1.25);
    hold on;
end
hold off;
title('Rough Heston smiles');
xlabel('Log-moneyness');
ylabel('Implied volatility');
hleg = legend();
title(hleg,'Expiry');

%% Compute option prices:
% Here we illustrate how one can also compute e.g. put option prices
% directly.

k = (-0.3:0.01:0.2)';
T = (0.2:0.2:3)';
[p,k_mat,T_mat] = model.GetPrices(k,T,true,'priceType','price',...
                                           'optionType','put');

p

% Also, if you want to use strikes instead of log-moneyness you need to compute
% them yourself like this:
F = model.s0.*exp((model.y.Eval(T) - model.q.Eval(T)).*T);
K_mat = exp(k_mat).*F.';

% Now we plot the result:
figure;
surf(K_mat,T_mat,p);
title('Put option prices under rough Bergomi');
xlabel('Strike');
ylabel('Expiration');
zlabel('Price');

%% A quick demo of some other functionality:
% Here we illustrate the various types of outputs that can be produced:
k = (-0.3:0.05:0.2)';
T = [0.5;0.75;1;2;3];
p = model.GetPrices(k,T,true,'priceType','price','optionType','put')
p = model.GetPrices(k,T,true,'priceType','price','optionType','call')
p = model.GetPrices(k,T,true,'priceType','implied_volatility')
p = model.GetPrices(k,T,true,'priceType','implied_volatility_surface')

% You can also request only out-of-the-money options (moneyness is here
% defined relative to the forward price):
[p,~,~,idxCall] = model.GetPrices(k,T,true,'priceType','price','optionType','otm');
p
idxCall

% The input vectors can also be interpreted row-by-row like this:
k = [0;-0.2;0;0.1;0.15];
T = [0.1;0.1;0.2;0.2;0.5];
cartProd = false;
p = model.GetPrices(k,T,cartProd,'priceType','implied_volatility')

%% Using more general curve objects:
% To define a more general forward variance curve we do as follows:
t_pts = [0.1;0.25;0.75;1;2];
xi_vals = [0.1;0.15;0.25;0.3;0.31].^2;
xi_curve = CurveClass('gridpoints',t_pts,'values',xi_vals);

% The curve looks like this:
t_eval = (0:0.01:3)';
xi_eval = xi_curve.Eval(t_eval);
figure;
plot(t_eval,xi_eval);
xlabel('Expiry');ylabel('Forward variance');title('Forward variance curve');
ylim([0,max(xi_eval)*1.1]);

% This can then be set in the object like this:
model.xi = xi_curve;

% Or we can simply define a new model as:
model = rHestonClass('H',0.1,'nu',0.3,'rho',-0.65,'xi',xi_curve,'s0',100);

% The same approach goes for the yield curve (stored in model.y) and the
% dividend yield curve (stored in model.q).

% The default interpolation method is 'flat'. You can also change this to
% 'linear' if you wish. An example (for a yield curve):
y_vals = [0.01;0.02;0.023;0.025;0.0255];
yield_curve = CurveClass('gridpoints',t_pts,'values',y_vals,...
                         'interpolation','linear');

y_eval = yield_curve.Eval(t_eval);
figure;
plot(t_eval,y_eval);
xlabel('Expiry');ylabel('Yield');title('Yield curve');
ylim([0,max(y_eval)*1.1]);

% As can be seen, the extrapolation is however still flat.

%% Settings:
% It is possible to change some of the settings of the pricing algorithm. More
% information can be found by reading the description of the rHestonClass
% and FourierPricerSettingsClass.
model.pricerSettings

