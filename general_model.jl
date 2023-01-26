using QuadGK
using Roots
using Plots
using Optim
using Distributions

# The function for tip given an expected tip (c) and number of honest bidders (n) found in equation X.
function t(F, c, θ, θ_threshold, n)
    if (θ < θ_threshold)
        return 0
    else
        integral_term, _ = quadgk(x -> F(x)^(n-1), 0, θ, rtol = 1e-8)
        return (1 / 2) * (integral_term - (n - 1) * c)
    end
end

function get_tip_function(F, c, θ_threshold, n)
    f(θ) = t(F, c, θ, θ_threshold, n)
    return f
end

function expected_tip(F, n)
    S(θ) = quadgk(x -> F(x)^(n-1), 0, θ, rtol = 1e-8)[1]

    #6.30
    u(x) = F(x)^(n-1) - F(x)^n
    six_thirty(x) = quadgk(u, x, 1, rtol =1e-8)[1] - 2 * S(x) / (n-1)

    θ_threshold = find_zero(six_thirty, (0, 1), Bisection())

    #6.24
    c = S(θ_threshold) / (n - 1)

    tip_function = get_tip_function(F, c, θ_threshold, n)

    return tip_function, c, θ_threshold
end

# Incomplete Beta Function
function b(x, α, β)
    u(t, α, β) = t^(α - 1) * (1 - t)^(β - 1)
    return quadgk(t -> u(t, α, β), 0, x, rtol =1e-8)[1]
end

# Beta Function
function B(α, β)
    u(t, α, β) = t^(α - 1) * (1 - t)^(β - 1)
    return quadgk(t -> u(t, α, β), 0, 1, rtol = 1e-8)[1]
end

# Beta distribution PDF
function f(x, α, β)
    return x^(α - 1) * (1 - x)^(β - 1) / B(α, β)
end


# Need to make sure this is at least 2
X = 2:100

plot(
    X, 
    [
        map(n -> n * expected_tip(x -> cdf(Beta(.5, .5), x), n)[2], X), 
        map(n -> n * expected_tip(x -> cdf(Beta(5, 1), x), n)[2], X),
        map(n -> n * expected_tip(x -> cdf(Beta(1, 3), x), n)[2], X),
        map(n -> n * expected_tip(x -> cdf(Beta(2, 2), x), n)[2], X),
        map(n -> n * expected_tip(x -> cdf(Beta(2, 5), x), n)[2], X),
    ],
    lw = 3, 
    labels = ["α = .5, β = .5" "α = 5, β = 1" "α = 1, β = 3" "α = 2, β = 2" "α = 2, β = 5"],
    title = "Expected Total Tip"
)

savefig("tip_plot_general")


X = 1:100
X = X / 100

plot(
    X, 
    [
        map(x -> f(x, .5, .5), X), 
        map(x -> f(x, 5, 1), X), 
        map(x -> f(x, 1, 3), X), 
        map(x -> f(x, 2, 2), X), 
        map(x -> f(x, 2, 5), X)], 
    lw = 3, 
    labels = ["α = .5, β = .5" "α = 5, β = 1" "α = 1, β = 3" "α = 2, β = 2" "α = 2, β = 5"], 
    title = "Beta Distribution PDFs"
)

savefig("beta_functions_general")

