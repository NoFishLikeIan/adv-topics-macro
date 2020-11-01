struct DetGrowthModel
    β::Float64
    α::Float64
    z::Float64
    ε::Float64
    k_size::Int
end

function equil_k(model::DetGrowthModel)
    model.z*(model.β / (1 + model.β))^(1 / (1 - model.α))
end