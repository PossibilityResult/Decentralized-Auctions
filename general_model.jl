using QuadGK
using Roots
using Plots
using Optim

# The function for tip given an expected tip (c) and number of honest bidders (n) found in equation X.
function t(F, c, θ, θ_threshold, n)
    F_to_n_minus_1(x) = F(x)^(n-1) 
    
    if (θ < θ_threshold)
        return 0
    else
        integral_term, _ = quadgk( F_to_n_minus_1, 0, θ, rtol = 1e-8)
        return (1 / 2) * ( integral_term - (n - 1) * c )
    end
end

function t̂(F, c, θ_threshold, n)
    f(θ) = t(F, c, θ, θ_threshold, n)
    return f
end

function param_finder(F, n)
    F_to_n_minus_1(x) = F(x)^(n-1) 
    F_to_n(x) = F(x)^n
    S(θ) = quadgk(F_to_n_minus_1, 0, θ, rtol = 1e-8)[1]

    #6.30
    function_to_integrate(x) = F_to_n_minus_1(x) - F_to_n(x)
    six_thirty(x) = quadgk(function_to_integrate, x, 1, rtol =1e-8)[1] - 2 * S(x) / (n-1)
    θ_threshold = find_zero(six_thirty, (0, 1), Bisection())

    #6.24
    c = S(θ_threshold)/(n-1)

    tip_function = t̂(F, c, θ_threshold, n)

    return tip_function, c, θ_threshold
end


function u(t, α, β)
    return t^(α - 1) * (1 - t)^(β - 1)
end

# Beta Function
function B(α, β)
    return quadgk(t -> u(t, α, β), 0, 1, rtol =1e-8)[1]
end

# Beta distribution PDF
function f(x, α, β)
    return x^(α - 1) * (1 - x)^(β - 1) / B(α, β)
end


expected_tip(n) = param_finder(F,n)[2]

# Need to make sure this is at least 2
X = 2:100


plot(X, X .* map(expected_tip, X), lw = 3, labels = ["nE[t]"], title = "Expected Total Tip")
savefig("tip_plot_general")
X = 1:100
X = X ./ 100

plot(
    [X, X, X, X, X], 
    [map(x -> F(x, .5, .5), X), map(x -> F(x, 5, 1), X), map(x -> F(x, 1, 3), X), map(x -> F(x, 2, 2), X), map(x -> F(x, 2, 5), X)], 
    lw = [3, 3, 3, 3, 3], 
    labels = ["α = .5, β = .5" "α = 5, β = 1" "α = 1, β = 3" "α = 2, β = 2" "α = 2, β = 5"], 
    title = "Beta Distribution PDFs"
)

savefig("beta_functions_general")

