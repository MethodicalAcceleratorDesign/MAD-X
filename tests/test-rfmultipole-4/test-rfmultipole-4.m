% Reads the file
M = load 'sectormap.dat';

% matrix J
J = [ 0 1 0 0  0 0 ; -1 0 0 0 0 0 ; 0 0 0 1 0 0 ; 0 0 -1 0 0 0 ; 0 ...
      0 0 0 0 1 ; 0 0 0 0 -1 0 ]

% extract the kick
K = M(1:6);

% extract the R matrix
R = reshape(M(6+(1:36)), 6, 6)

% extract the T matrix
T = {};
I = zeros(1, 216);
index=1;
for k=0:5
    for j=0:5
        for i=0:5
            I(index++) = 6 + 36 + 1 + (i*36) + (j*6) + k;  
        end
    end
end

for k=0:5
    T{k+1} = zeros(6,6);
    for i=0:5
        for j=0:5
            T{k+1}(i+1,j+1) = M(I((i*36) + (j*6) + k + 1)); 
        end
    end
end

disp('Test 1: transpose(R) * J * R == J')

R' * J * R == J

disp(['Test 2: transpose(Tk) * J * R + transpose(R) * J * Tk == 0 ' ...
      'for k=1..6']);

for k=1:6
    k
    T{k}' * J * R + R' * J * T{k}
end

keyboard
