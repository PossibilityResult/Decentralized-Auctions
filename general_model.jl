using QuadGK
using Roots
using Plots
using Optim



# The function for tip given an expected tip (c) and number of honest bidders (n) found in equation X.
function t(F, c, θ, n)
    F_to_n(x) = F(x)^n 
    F_to_n_minus_1(x) = F(x)^(n-1) 
    S(θ_threshold), _ = quadgk(F_to_n_minus_1, 0, θ_threshold, rtol = 1e-8)
    g(θ_threshold) = S(1) - quadgk(F_to_n, θ_threshold, 1, rtol = 1e-8)[1] - (2n - (1 - θ_threshold)) * S(θ_threshold) / (n - 1)
    
    θ_threshold = find_zero(g, (0, 1), Bisection())

    if (θ < θ_threshold)
        return 0
    else
        integral_term, _ = quadgk( F_to_n_minus_1, 0, θ, rtol = 1e-8)
        return (1 / 2) * ( integral_term - (n - 1) * c )
    end
end

function t̂(F, c, n)
    f(θ) = t(F, c, θ, n)
    return f
end

# The function for expected tip given the number of honest bidders (n) found in equation X.
function expected_tip(F, n)
    function g(c)
        I, _ = quadgk(t̂(F, c, n), 0, 1, rtol = 1e-8)
        return I - c
    end
    return find_zero(g, (0, 1), Bisection())
end

function expected_tip̂(n)
    F(x) = x
    return expected_tip(F, n)
end



# --- PLOTTING THE EXPECTED TOTAL AMOUNT TIPPED AND UPPER BOUND ---



function g(θ_threshold)
    F(x) = x
    n = 8
    F_to_n(x) = F(x)^n
    F_to_n_minus_1(x) = F(x)^(n-1) 
    S(θ_threshold), _ = quadgk(F_to_n_minus_1, 0, θ_threshold, rtol = 1e-8)
    return S(1) - quadgk(F_to_n, θ_threshold, 1, rtol = 1e-8)[1] - (2n - (1 - θ_threshold)) * S(θ_threshold) / (n - 1)
end

# Column vector of the number of honest bidders to iterate through and plot.
X = 1:150
X = X ./ 100

# Plotting the total expected tip and an upper bound of it.
plot(g, X, label = "g(θ_threshold)")
#plot(X, [X .* map(expected_tip̂, X)], lw = [3, 3], labels = ["nE[t]" "1 / ((n-1) * sqrt(n))"], title = "Expected Total Tip")