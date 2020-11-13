using Dolo, Distances


model = yaml_import("src/week-three/models/sgm.yaml")

no_shocks = Dict(:z => [model.calibration.flat[:z]])

T = 101 # Simulation periods
p_T = 10 # Plotting periods
shocks = Dict(:z => [fill(0.2, 2); 0.4])


sol_ref = perfect_foresight(
    model, shocks, 
    T=T, complementarities=false)
