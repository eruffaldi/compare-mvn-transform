function r = symmtx(name,size,cond)

r = sym(name,size);
if isempty(cond) == 0
    try
        r = sym(r,cond);
    catch me
        assume(r,cond);
    end
end
