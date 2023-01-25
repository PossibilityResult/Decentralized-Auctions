using QuadGK
using Roots
using Plots

# Point at which expected utility of tipping is greater than 0.
# Before this point, bidders tip 0.
θ_threshold(c, n) = ((n^2 - n) * c)^(1 / n)


# The function for tip given an expected tip (c) and number of honest bidders (n) found in equation X.
function t(c, θ, n)
    if (θ < θ_threshold(c, n))
        return 0
    else
        return (1 / 2)*(θ^n / n - (n - 1) * c)
    end
end

function t̂(c, n)
    f(θ) = t(c, θ, n)
    return f
end

# The function for expected tip given the number of honest bidders (n) found in equation X.
function expected_tip(n)
    function g(c)
        I, est = quadgk(t̂(c, n), 0, 1, rtol = 1e-8)
        return I - c
    end
    return find_zero(g, (0, 1), Bisection())
end


function g(c)
    I, est = quadgk(t̂(c, n), 0, 1, rtol = 1e-8)
    return I - c
end

print(find_zero(g, (0, 1), Bisection()))

# --- PLOTTING THE EXPECTED TOTAL AMOUNT TIPPED AND UPPER BOUND ---

# Column vector of the number of honest bidders to iterate through and plot.
X = 1:100

# The upper bound on the expected total amount tipped.
s(n) = 1 / ((n - 1) * sqrt(n))

# Plotting the total expected tip and an upper bound of it.
plot(X, [X .* map(expected_tip, X), map(s, X)], lw = [3, 3], labels = ["nE[t]" "1 / ((n-1) * sqrt(n))"], title = "Expected Total Tip and Upper Bound")
