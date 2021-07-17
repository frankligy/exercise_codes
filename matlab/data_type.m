%{ basic data types }%

a = 1;
b = 'hello';
c = "hello";
d = [1 2 3];  % row vector
e = [1;2;3];  % column vector
f = [1 2 3; 4 5 6];
g = ["hello" "yes"];
h = {'hello' 34 [1 2 3];'yes' 45 [4 5 66]}; % cell array
i = struct('name','frank','gpa',3.9,'grade',3);
j(1,2,2) = i;
j(1,1,2) = i;

a + 1

%{ how to access element in array }%

g(1,2)
g(1,end)  % end means the last index in the axis
h{1,2}
g(:,1:1:2) % demo slicing start:step:end
g(:,2:-1:1) % slicing in reversed orientation