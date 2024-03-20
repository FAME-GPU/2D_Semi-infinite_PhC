function is_parallel = isParallel(x, y)

x = x / norm(x);
y = y / norm(y);

if norm(x - y, 1) < 1e-14
    is_parallel = 1;
else
    is_parallel = 0;
end

end