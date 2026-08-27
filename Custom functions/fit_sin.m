function [fitobj,gof] = fit_sin(t,y, starting_points)

v = ver;

fun = @(x)x(1)*sin(x(2)*t + x(3)) + x(4) - y;
options = optimoptions(@lsqnonlin,'Display','off');

if str2double(v(end).Date(end-3:end))>=2024
    x = lsqnonlin(fun,starting_points,[],[],[1 0 0 -1],0,[],[],[],options);
else
    x = lsqnonlin(fun,starting_points,[],[],options);
end
predicted_value = x(1)*(sin(x(2)*t + x(3))) + x(4);
sserr       =sum((predicted_value-y).^2);
sstot       = sum((y - mean(y)).^2);
gof         = [];
gof.sse     = sserr;
gof.sst     = sstot;
gof.rmse    = sqrt(sserr/numel(y));
gof.rsquare = 1 - sserr/sstot;
x           = num2cell(x);
mdl = fittype('a*sin(b*x + c) + d','indep','x');
fitobj = cfit(mdl, x{:});
end