using JuMP
using Ipopt
using LinearAlgebra
using Plots
using DrWatson
using Dates

time_now = Dates.format(now(), "Y-mm-dd-HH:MM:SS")
outputdir(args...) = projectdir("output", args...)
outputdir_subpath = outputdir(time_now)
mkdir(outputdir_subpath)

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_cpu_time", 60.0)
set_optimizer_attribute(model, "print_level", 0)


@variable(model, Z[1:31] ≥ 0)

D = zeros(30,31)
D[diagind(D)] .= -1
D[diagind(D,1)] .= 1
f = sum((D * Z).^2) * 30 / sqrt(1.6^2 -1)
@objective(model, Min, f)

A = zeros(30,31)
A[diagind(D)] .= 1
A[diagind(D,1)] .= 1

# @constraints(model, begin
#     Z[1]  == 0
#     Z[31] == 0 
# end)

# @constraints(model, begin
#     Z[1]  == 0
#     Z[31] == 0 
#     Z[16] ≥  0.21
# end)

# @constraints(model, begin
#     Z[1]  == 0
#     Z[31] == 0
#     sum(A * Z) == 15/2
# end)

@constraints(model, begin
    Z[1]  == 0
    Z[31] == 0
    Z[16] ≥  0.21
    sum(A * Z) == 15/2
end)


print(model)
optimize!(model)
println(solution_summary(model))
# @show termination_status(model)
# @show primal_status(model)
# @show dual_status(model)
# @show objective_value(model)
# @show value.(Z)

x = range(-1,1,length=31)
output = value.(Z)

plot([x, x],[output,-output], xlabel="x", ylabel="z", legend=false)

savefig(outputdir(outputdir_subpath, "foil.png"))