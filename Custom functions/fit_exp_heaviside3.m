function fit_hs_exp = fit_exp_heaviside3(t,y, starting_points)


fited_param = lsqnonlin(@myfun,starting_points);

fit_hs_exp.fited_param = fited_param;
x = fited_param;
fited_hs = heaviside(t-x(1))*(x(3)-x(2)) +x(2);
predicted_value = (x(4)*exp(-t.*(fited_hs*(x(3)-x(2)) +x(2)))  + x(5) + fited_hs*(x(4)*(exp(-x(1)*x(2))-exp(-x(1)*x(3)))))';
fit_hs_exp.predicted_value = predicted_value;


    function F = myfun(x)
        hs_function = heaviside(t-x(1));
        F = (x(4)*exp(-t.*(hs_function*(x(3)-x(2)) +x(2)))  + x(5) + hs_function*(x(4)*(exp(-x(1)*x(2))-exp(-x(1)*x(3)))) + - y)';
    end

end