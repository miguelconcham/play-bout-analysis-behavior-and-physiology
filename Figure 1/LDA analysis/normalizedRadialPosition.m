function pos = normalizedRadialPosition(x, y, x_inner, y_inner, x_outer, y_outer, center)

    % Compute direction vector from center to point

    if  inpolygon(x,y,x_inner,x_inner)
        pos = 0
    elseif ~inpolygon(x,y,x_outer,y_outer)
        pos = 1;
    else
    v = [x - center(1), y - center(2)];
    v_unit = v / norm(v);

    % Ray from center in direction v
    ray_length = 1e3; % Large enough
    ray_end = center + ray_length * v_unit;

    % Intersect ray with inner and outer shapes
    [xi_in, yi_in] = polyxpoly([center(1), ray_end(1)], [center(2), ray_end(2)], x_inner, y_inner);
    [xi_out, yi_out] = polyxpoly([center(1), ray_end(1)], [center(2), ray_end(2)], x_outer, y_outer);

    % Choose closest intersections along ray direction
    d_inner = min(vecnorm([xi_in - center(1), yi_in - center(2)], 2, 2));
    d_outer = min(vecnorm([xi_out - center(1), yi_out - center(2)], 2, 2));

    % Distance of point from center
    d_point = norm(v);

    % Normalize
    pos = (d_point - d_inner) / (d_outer - d_inner);
    pos = min(max(pos, 0), 1); % Clamp to [0,1]
    end
end