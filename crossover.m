function child = crossover(parent1, parent2)

n = size(parent1,2);
child = [];
k = floor(0.3*n);
t1 = parent1;
t2 = parent2;

while n > 0
    t11_size = min(k, size(t1,2));
    t11_start = randi(n - t11_size + 1);
    child = [child t1(1,t11_start:t11_start+t11_size-1)];
    for i = 1:size(child,2)
        t1 = t1(t1 ~= child(1,i));
        t2 = t2(t2 ~= child(1,i));
    end
    x = t1;
    t1 = t2;
    t2 = x;
    n = size(t1,2);
end
