
include(joinpath(dirname(Base.active_project()), "src", "utils", "interval.jl"))

using .Intervals

using JCheck
import JCheck: generate, shrinkable, shrink

# Generators

function sample(rng, items, weights) 
    @assert length(items) == length(weights)
    threshold = rand(rng) * sum(weights)
    i = findfirst(cumsum(weights) .>= threshold)
    return items[i]
end

function generate_bound(rng)
    return sample(
        rng,
        [(rand(rng) - 0.5) * 100., -Inf, Inf, -0., 0.],
        [10, 1, 1, 1, 1])
end

function generate_interval(rng)
    candidate = Interval(generate_bound(rng), generate_bound(rng))
    if Intervals.is_valid(candidate) 
        return candidate
    else
        return generate_interval(rng)
    end
end

function generate(rng, ::Type{Interval}, n :: Int)
    return [generate_interval(rng) for _ in 1:n]    
end


# Tests

tests = Quickcheck("Interval Tests")

@add_predicate tests "Addition is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(a + b))

@add_predicate tests "Substraction is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(a - b))

@add_predicate tests "Negation is valid" ((a :: Interval) -> 
    Intervals.is_valid(-a))

@add_predicate tests "Multiplication is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(a * b))
    
@add_predicate tests "Inversion is valid" ((a :: Interval) -> 
    Intervals.is_valid(inv(a)))

@add_predicate tests "Division is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(a / b))

@add_predicate tests "Minimum is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(min(a, b)))

@add_predicate tests "Maximum is valid" ((a :: Interval, b :: Interval) -> 
    Intervals.is_valid(max(a, b)))

@add_predicate tests "Sine is valid" ((a :: Interval) -> 
    Intervals.is_valid(sin(a)))

@add_predicate tests "Cosine is valid" ((a :: Interval) -> 
    Intervals.is_valid(cos(a)))

@add_predicate tests "Exponential is valid" ((a :: Interval) -> 
    Intervals.is_valid(exp(a)))

@add_predicate tests "Square root is valid" ((a :: Interval) -> 
    Intervals.is_valid(sqrt(a)))

@quickcheck tests


