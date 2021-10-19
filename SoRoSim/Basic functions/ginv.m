function g = ginv(g)
g          = [g(1:3,1:3)' -g(1:3,1:3)'*g(1:3,4);0 0 0 1];
end

