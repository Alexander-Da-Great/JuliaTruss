# insert statements to call in functions from the helper files in this folder
using MAT
using Plots
using TickTock
include("static_solver.jl")
# step 1: read in the data
threshold = 0.1;
tick()
#=
% ==============================================
% CHANGE THIS VARIABLE TO TEST DIFFERENT TRUSSES
% truss = 1 for medium truss      
% truss = 2 for large truss
% truss = 3 for small truss with thermal effects
% truss = 4 for large truss with thermal effects
truss = 4;
% ==============================================
/Users/alexandergreer/Documents/Julia/FEM/
=#
truss = 4
start = pwd()
if truss == 1
    #cd("data")
    vars = matread(start*"/data/med_truss.mat")
    sols = matread(start*"/data/med_truss_answers.mat")
    min_u = 4;
    min_F = 4;
    min_sigma = 0;
elseif truss == 2
    vars = matread(start*"/data/large_truss.mat")
    sols = matread(start*"/data/large_truss_answers.mat")
    min_u = 12;
    min_F = 21*3;
    min_sigma = 0;
elseif truss == 3
    vars = matread(start*"/data/small_truss_thermal.mat")
    sols = matread(start*"/data/small_truss__thermal_answers.mat")
    min_u = 6;
    min_F = 2;
    min_sigma = 0;
elseif truss == 4
    vars = matread(start*"/data/large_truss_thermal.mat")
    sols = matread(start*"/data/large_truss_thermal_answers.mat")
    min_u = 12;
    min_F = 21*3;
    min_sigma = 0;
end
nodetable = vars["nodetable"]
connecttable = vars["connecttable"] 
dof = convert(Int64,vars["dof"])
min_score = min_u + min_F + min_sigma;

#sol = test_print(nodetable, connecttable, dof)

sol = static_solver(nodetable, connecttable, dof)
#print(sol[1])
u = sol[1]; # displacements
plotlyjs()
# now we plot
plt = plot()
for m=1:size(connecttable,1)
    plot!([nodetable[Int64(connecttable[m,2]),2];nodetable[Int64(connecttable[m,3]),2]],
          [nodetable[Int64(connecttable[m,2]),3];nodetable[Int64(connecttable[m,3]),3]],
          [nodetable[Int64(connecttable[m,2]),4];nodetable[Int64(connecttable[m,3]),4]], legend = false, label = ["initial"], linecolor = :black);
    point1 = u[Int64(connecttable[m,2]),:]
    point2 = u[Int64(connecttable[m,3]),:]
    #print(nodetable[Int64(connecttable[m,2]),2] - point1[1])
    plot!([nodetable[Int64(connecttable[m,2]),2] + point1[1];nodetable[Int64(connecttable[m,3]),2] + point2[1]],
          [nodetable[Int64(connecttable[m,2]),3] + point1[2];nodetable[Int64(connecttable[m,3]),3] + point2[2]],
          [nodetable[Int64(connecttable[m,2]),4] + point1[3];nodetable[Int64(connecttable[m,3]),4] + point2[3]], legend = false, label = ["final"], linecolor = :red);

    # insert code to include the stresses and forces in the members

end
plot!(title = "Thermo-structural Truss Sim",legend = :none, label = ["initial" "final"], aspect_ratio = 1)
display(plt)
#png(plt,"truss")
tock()
