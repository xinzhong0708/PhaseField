load temp
X  = X_Fit(100,:);
dc = 1e-5;
X1 = X;  X1(4)=X1(4)+dc;
% Forward pass
Z1 = tanh(X*W1' + b1');       % N x hidden1
Z2 = tanh(Z1*W2' + b2');      % N x hidden2
f  = Z2*W3' + b3';            % N x 1
% Forward pass
Z1 = tanh(X1*W1' + b1');       % N x hidden1
Z2 = tanh(Z1*W2' + b2');      % N x hidden2
f1 = Z2*W3' + b3';            % N x 1
(f1-f)/dc

% For efficiency, compute with matrix multiplications
mu = zeros(size(X));  % N x d
for i = 1:size(X,1)
    J1 = diag(1 - Z1(i,:).^2) * W1;   % hidden1 x d
    J2 = diag(1 - Z2(i,:).^2) * W2;   % hidden2 x hidden1
    mu(i,:) = (W3 * J2 * J1)';         % 1 x d
end
mu


% Get the number of input dimensions (d)
d = size(X, 2);
Hessian = zeros(d, d);

% Compute the baseline gradient for the original input X
mu_base = compute_mu(X, W1, b1, W2, b2, W3, b3);

% Loop through each input dimension to compute the Hessian columns
for j = 1:d
    % Create a perturbed input vector X_perturbed
    X_perturbed = X;
    X_perturbed(j) = X_perturbed(j) + dc;

    % Compute the gradient for the perturbed input
    mu_perturbed = compute_mu(X_perturbed, W1, b1, W2, b2, W3, b3);

    % Compute the j-th column of the Hessian using finite differences
    Hessian(:, j) = (mu_perturbed - mu_base)' / dc;
end

disp('Hessian matrix:');
disp(Hessian);



% Function to compute the gradient (your mu variable)
function mu_val = compute_mu(X_input, W1, b1, W2, b2, W3, b3)
    N = size(X_input, 1);
    mu_val = zeros(N, size(X_input, 2));

    for i = 1:N
        Z1_i = tanh(X_input(i,:)*W1' + b1');
        Z2_i = tanh(Z1_i*W2' + b2');
        
        J1 = diag(1 - Z1_i.^2) * W1;   % hidden1 x d
        J2 = diag(1 - Z2_i.^2) * W2;   % hidden2 x hidden1
        mu_val(i,:) = (W3 * J2 * J1)'; % 1 x d
    end
end