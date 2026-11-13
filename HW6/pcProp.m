function [t_out, states_out, poincare] = pcProp(x0, u0, v0, mu, max_crossings)
    % Propagate trajectory and collect Poincare sections
    
    y0 = [x0; 0; 0; u0; v0; 0];
    
    % Integration parameters
    tspan = [0 1000]; % Max time
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, ...
                     'Events', @(t,y) xAxisCross(t,y));
    
    poincare = [];
    t_all = [];
    states_all = [];
    y_current = y0;
    t_current = 0;
    
    for i = 1:max_crossings
        [t, y, te, ye, ie] = ode113(@(t,y)EarthMoonCR3BP(t,y,mu), ...
                                     [t_current, t_current+100], y_current, options);
        
        if isempty(te)
            break;
        end
        
        t_all = [t_all; t];
        states_all = [states_all; y];
        
        % Record crossing
        poincare = [poincare; ye(end,1), ye(end,3), ye(end,4), te(end)];
        
        % Continue from crossing
        y_current = ye(end,:)';
        t_current = te(end);
    end
    
    t_out = t_all;
    states_out = states_all;
end