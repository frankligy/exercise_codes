% if statement
a = 5;
if a < 30
    disp('small')
elseif a > 30
    disp('big')
else
    disp('equal')
end

% for loop
for i = 1:a
    disp(i)
end

% switch statement
switch a
    case 1
        disp('it is 1')
    case 2
        disp('it is 2')
    case 3
        disp('it is 3')
    case 4
        disp('it is 4')
    case 5
        disp('correct')
end

% while statement
i = 1;
while i <= 5
    disp(i)
    i = i + 1;
end

% custom function
return_value = function_name(1,2);
function return_value = function_name(arg1,arg2)
    return_value = arg1 + arg2;
end