function x = analytical_solution(t)
    n = 0;
    t1 = 2*(pi*n + atan(1/15*(-1 - 22*(2/(101 + 15*sqrt(69)))^(1/3) + 2^(2/3)*(101 + 15*sqrt(69))^(1/3))));
    x11 = 0.5*sin(t1); x21 = 0.5*cos(t1);
end
