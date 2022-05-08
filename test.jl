# This is a test script to aid in the transliation of code from MATLAB to julia
using LinearAlgebra
using MAT
A = rand(3,3)
norm(A)
start = pwd()
print(start)
A'

F = [0 0
50 -86
  NaN NaN
  NaN NaN];

A = rand(5, 5)
b = A .> 0.5
v = [1 1 1 0 0]
print(A[v'*v])

function test_print(nodetable, connecttable, dof)
    println(nodetable)
    print("1")
    println(connecttable)
    print("2")
    println(dof)
    print("3")
    return [1,2,3]
end