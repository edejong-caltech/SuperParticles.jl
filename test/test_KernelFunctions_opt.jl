using SuperParticles.KernelFunctions
using JET: @test_opt

rtol = 1e-3

# Constant kernel
# Initialization
@test_opt ConstantKernelFunction(FT(π))

# Evaluations
x = 2.0
y = 4.0
kernel = ConstantKernelFunction(FT(π))
@test_opt kernel(x, y)


# Linear kernel
# Initialization
@test_opt LinearKernelFunction(FT(π))

# Evaluations
x = 2.0
y = 4.0
kernel = LinearKernelFunction(FT(π))
@test_opt kernel(x, y)