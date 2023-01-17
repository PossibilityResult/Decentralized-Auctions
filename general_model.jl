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
    function_to_integrate(x) = F_to_n_minus_1(x)-F_to_n(x)
    six_thirty(x) = quadgk(function_to_integrate,x,1,rtol =1e-8)[1] - 2*S(x)/(n-1)
    θ_threshold = find_zero(six_thirty, (0, 1), Bisection())

    #6.24
    c = S(θ_threshold)/(n-1)

    tip_function = t̂(F, c, θ_threshold, n)

    return tip_function, c, θ_threshold
end

function F(x)
    if x<.5
        return(2*x^2)
    else
        return(2*x^2 - (2x-1)^2)
    end
end

n = 10

expected_tip(n) = param_finder(F,n)[2]

#need to make sure this is at least 2
X = 2:100


plot(X, X .* map(expected_tip, X), lw = [3, 3], labels = ["nE[t]"], title = "Expected Total Tip")
savefig("tip_plot_general")



# F(x) = x 
# n = 10
# c = 0.002259854780578178
# θ_threshold  =  0.8527707666784111

# plot(X, map(S,X), fmt = :pdf)
# savefig("S_plot")

#6.30

# function six_thirty(θ_threshold)
#     function_to_integrate(x) = F_to_n_minus_1(x)-F_to_n(x)
#     left_hand_side = quadgk(function_to_integrate,θ_threshold,1,rtol =1e-8)[1]
#     return(left_hand_side - 2*S(θ_threshold)/(n-1))
# end

# print(find_zero(six_thirty, (0, 1), Bisection()))

# plot(X, map(six_thirty,X), fmt = :pdf)
# savefig("6_30_plot")

# # The function for expected tip given the number of honest bidders (n) found in equation X.
# function expected_tip(F, n)
#     function g(c)
#         I, _ = quadgk(t̂(F, c, n), 0, 1, rtol = 1e-8)
#         return I - c
#     end
#     return find_zero(g, (0, 1), Bisection())
# end

# function expected_tip̂(n)
#     F(x) = x
#     return expected_tip(F, n)
# end



# --- PLOTTING THE EXPECTED TOTAL AMOUNT TIPPED AND UPPER BOUND ---



# function g(θ_threshold)
#     F(x) = x
#     n = 8
#     F_to_n(x) = F(x)^n
#     F_to_n_minus_1(x) = F(x)^(n-1) 
#     S(θ_threshold), _ = quadgk(F_to_n_minus_1, 0, θ_threshold, rtol = 1e-8)
#     return S(1) - quadgk(F_to_n, θ_threshold, 1, rtol = 1e-8)[1] - (2n - (1 - F(θ_threshold)) ) * S(θ_threshold) / (n - 1)
# end

# function h(θ_threshold)
#     F(x) = x
#     n = 8
#     F_to_n(x) = F(x)^n
#     F_to_n_minus_1(x) = F(x)^(n-1) 
#     S(θ_threshold), _ = quadgk(F_to_n_minus_1, 0, θ_threshold, rtol = 1e-8)
#     return S(1) - S(θ_threshold) - quadgk(F_to_n, θ_threshold, 1, rtol = 1e-8)[1] - (2 + (n - 1) * (1 - F(θ_threshold)) ) * S(θ_threshold) / (n - 1)
# end

# function q(θ_threshold)
#     F(x) = x
#     n = 8
#     F_to_n(x) = F(x)^n
#     F_to_n_minus_1(x) = F(x)^(n-1) 
#     S(θ_threshold), _ = quadgk(F_to_n_minus_1, 0, θ_threshold, rtol = 1e-8)
#     return S(1) - S(θ_threshold) - quadgk(F_to_n, θ_threshold, 1, rtol = 1e-8)[1] - (1 + n) * S(θ_threshold) / (n - 1)
# end




# Column vector of the number of honest bidders to iterate through and plot.
# X = 1:100
# X = X ./ 100

# # Plotting the total expected tip and an upper bound of it.
# plot(X, [g, h, q], labels = ["g(θ_threshold" "h(θ_threshold)" "q(θ_threshold)"])

#n = 8
#F(x) = x
#F_to_n_minus_1(x) = F(x)^(n-1) 
#S(x) = quadgk(F_to_n_minus_1, 0, x, rtol = 1e-8)[1]
#plot(S, X, label = "S(x)")

#plot(X, [X .* map(expected_tip̂, X)], lw = [3, 3], labels = ["nE[t]" "1 / ((n-1) * sqrt(n))"], title = "Expected Total Tip")
