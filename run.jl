using JuMP
using Ipopt
using LinearAlgebra
using Plots
using DrWatson
using Dates
using DataFrames, CSV

time_now = Dates.format(now(), "Y-mm-dd-HH:MM:SS")
outputdir(args...) = projectdir("output", args...)
outputdir_subpath = outputdir(time_now)
mkdir(outputdir_subpath)
obj_values = []

# Matrices for differencing and adding two consecutive elements of a vector
D = zeros(30,31)
D[diagind(D)] .= -1
D[diagind(D,1)] .= 1

A = zeros(30,31)
A[diagind(D)] .= 1
A[diagind(D,1)] .= 1

function add_constraint_case1(model)
    Z = model[:Z]
    @constraints(model, 
    begin
        Z[1]  == 0
        Z[31] == 0 
    end
    )
end

function add_constraint_case2(model)
    Z = model[:Z]
    @constraints(model, 
    begin
        Z[1]  == 0
        Z[31] == 0
        Z[16] ≥  0.21
    end
    )
end

function add_constraint_case3(model)
    Z = model[:Z]
    @constraints(model, 
        begin
            Z[1]  == 0
            Z[31] == 0
            sum(A * Z) == 15/2
        end
    )
end

function add_constraint_case4(model)
    Z = model[:Z]
    @constraints(model, 
    begin
        Z[1]  == 0
        Z[31] == 0
        Z[16] ≥  0.21
        sum(A * Z) == 15/2
    end
    )
end

constraint_fns = [add_constraint_case1, add_constraint_case2, add_constraint_case3, add_constraint_case4]
constraint_names = ["Closed Body", "Closed Body + Min Thickness", "Closed Body + Fixed Area", "All Constraints"]
x = range(-1,1,length=31)
plots = []
areas = []
for (i,(c_fn, c_name)) in enumerate(zip(constraint_fns, constraint_names))
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", 60.0)
    set_optimizer_attribute(model, "print_level", 0)

    @variable(model, Z[1:31] ≥ 0)
    @objective(model, Min,
        sum((D * Z).^2) * 30 / sqrt(1.6^2 -1)
    )

    c_fn(model)
    optimize!(model)

    output = value.(Z)
    push!(obj_values, objective_value(model))
    push!(areas, sum(A * output) * 2 / 30)

    push!(plots,
        plot([x, x],[output,-output],
            xlabel="x", ylabel="z", legend=false,
            ylims=[-0.2, 0.2],
            title="Case $i: $c_name", titlefontsize=8,
            grid = true, gridalpha = 0.5, 
        )
    )

    if c_name == "Closed Body"
        plot([x, x],[output,-output],
            xlabel="x", ylabel="z", legend=false,
            title="Case $i: $c_name", titlefontsize=8,
            grid = true, gridalpha = 0.5, 
        )
        savefig(outputdir(outputdir_subpath, "case1.png"))
    end
end

plot(plots...)
savefig(outputdir(outputdir_subpath, "all.png"))

df = DataFrame(Constraints=constraint_names, Cd=obj_values, Areas = areas)
fout = outputdir(outputdir_subpath, "data.csv")
CSV.write(fout, df)