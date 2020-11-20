function toupper(x::String, i::Int)::String
    x[1:i - 1] * uppercase(x[i:i]) * x[i + 1:end]
end