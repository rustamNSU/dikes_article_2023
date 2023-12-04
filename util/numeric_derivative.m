function [v, tv] = numeric_derivative(s, t)
    if not(length(s) == length(t))
        error('central_difference(s, t): different length of s and t');
    end
    n = length(s);
    v = zeros(1, n-1);
    tv = zeros(1, n-1);
    for i = 1:n-1
        v(i) = (s(i+1)-s(i))/ (t(i+1)-t(i));
        tv(i) = t(i+1);
    end
end

