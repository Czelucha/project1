%Define Eta function
function f = Eta(r)
%Function Eta evaluates r^2log(r) if r > 0, else it returns 0
    if r == 0
        f = 0;
    else
        f = r^2 * log(r);  
    end
end