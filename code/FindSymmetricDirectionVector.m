function symDir = FindSymmetricDirectionVector(v, a, b, c, n)
%Calculate the symmetric direction vector at a point about another direction vector,
%where the point is also constrainted to a plane

bSys = [
    dot(v,(b-c));
    dot(a,b) - sumsqr(a);
    dot(n,(a-b));
    ];

A = [
    v';
    c';
    n';
    ];

symVec = A\bSys;

symDir = symVec./norm(symVec);
end

