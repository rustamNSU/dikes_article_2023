function [v, tv] = central_difference(s, t)
    if not(length(s) == length(t))
        error('central_difference(s, t): different length of s and t');
    end
    n = length(s);
    v = zeros(1, n-4);
    tv = zeros(1, n-4);
    for i = 1:n-4
        v(i) = (s(i+4)-s(i))/ (t(i+4)-t(i));
        tv(i) = t(i+4);
    end
end

