var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#FourierFlows.jl-Documentation-1",
    "page": "Home",
    "title": "FourierFlows.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "FourierFlows provides solvers for partial differential equations on doubly-periodic domains using Fourier-based pseudospectral methods. A central intent of the software's design is also to provide a framework for writing new, fast solvers for new physical problems. The code is written in Julia."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "FourierFlows is a registered package and so can be installed via Pkg.add.Pkg.add(\"FourierFlows\")For now, this package supports Julia 0.6. Support for version 0.7 is on its way."
},

{
    "location": "index.html#Basic-Notation-1",
    "page": "Home",
    "title": "Basic Notation",
    "category": "section",
    "text": "The code solves partial differential equations of the general form:partial_t u = mathcalLu + mathcalN(u) The mathcalLu part is time-stepped forward using an implicit scheme; the mathcalN(u) part is time-stepped forward using an explicit scheme.The coefficients for the linear operator mathcalL are stored in array LC. The term mathcalN(u) is computed for by calling the function calcN."
},

{
    "location": "index.html#Source-code-organization-1",
    "page": "Home",
    "title": "Source code organization",
    "category": "section",
    "text": "The code is divided along conceptual lines into problem-agnostic and problem-specific components. Files that contain problem-agnostic parts of the code are stored in /src. Files in /src define the domain, 'AbstractTypes' that supertype problem-specific types, and time-stepper types and routines. Problem-specific modules are stores in /src/physics.Here's an overview of the code structure:/src/\nFourierFlows.jl\nDefines supertyping AbstractParams, AbstractGrid, etc.\nDefines a Problem type to organize the grid, vars, params,   equation, and timestepper into a single structure.\nIncludes all sources files and physics files.\ntimesteppers.jl: defines modules and stepforward! routines for   various time-steppers. Current implemented time-steppers are:\nForward Euler (+ Filtered Forward Euler)\n3rd-order Adams-Bashforth (AB3)\n4th-order Runge-Kutta (RK4)\n4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)\n(+ Filtered ETDRK4)\nphysics/\ntwodturb.jl: Defines a TwoDTurb module that provides a       solver for the two-dimensional vorticity equation.\nbarotropicqg.jl: Defines a BarotropicQG module that provides       several solvers for the barotropic QG model that permit beta,       topography, beta + topography, and forcing.\ntwomodeboussinesq.jl: Defines a TwoModeBoussinesq module       that provides solvers for a two-mode truncation of the       rotating, stratified Boussinesq equation.\ntraceradvdiff.jl: Defines a TracerAdvDiff module that       provides a solver for a two-dimensional and periodic tracer       field in a given 2D flow (u, w), which can be an arbitrary       function of x, z, and t.\ntracerpatcheqn.jl: ..."
},

{
    "location": "index.html#Writing-fast-solvers-1",
    "page": "Home",
    "title": "Writing fast solvers",
    "category": "section",
    "text": "The performance-intensive part of the code involves just two functions: the timestepping scheme stepforward!, and the function calcN! that calculates the nonlinear part of the given equation's right-hand side. Optimization of these two functions for a given problem will produce the fastest possible code."
},

{
    "location": "index.html#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "An example script that simulates decaying two-dimensional turbulence reproducing the results of the paper byMcWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. J. Fluid Mech., 146, 21-43.is found in examples/twodturb/McWilliams.jl."
},

{
    "location": "index.html#Future-work-1",
    "page": "Home",
    "title": "Future work",
    "category": "section",
    "text": "The code is in the chaos stage of development. A main goal for the future is to permit the use of shared memory parallelism in the compute-intensive routines (shared-memory parallelism provided already by FFTW/MKLFFT, but is not yet native to Julia for things like element-wise matrix multiplication, addition, and assignment). This feature may possibly be enabled by Intel Lab's ParallelAccelerator package."
},

{
    "location": "index.html#Authors-1",
    "page": "Home",
    "title": "Authors",
    "category": "section",
    "text": "FourierFlows is currently being developed by Gregory L. Wagner and Navid C. Constantinou."
},

{
    "location": "index.html#DocStrings-1",
    "page": "Home",
    "title": "DocStrings",
    "category": "section",
    "text": "Pages = [\n    \"man/docstrings.md\",\n    ]\nDepth = 2"
},

{
    "location": "index.html#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\n    \"man/docstrings.md\",\n    ]"
},

{
    "location": "man/docstrings.html#",
    "page": "Functions exported from FourierFlows:",
    "title": "Functions exported from FourierFlows:",
    "category": "page",
    "text": ""
},

{
    "location": "man/docstrings.html#Base.resize!-Tuple{FourierFlows.AbstractDiagnostic,Int64}",
    "page": "Functions exported from FourierFlows:",
    "title": "Base.resize!",
    "category": "Method",
    "text": "resize!(diag, newnum)\n\nResize the Diagnostic data and time arrays to length newnum. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.increment!-Tuple{AbstractArray}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.increment!",
    "category": "Method",
    "text": "increment!(diags)\n\nIncrement the array of Diagnostics diags. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.increment!-Tuple{FourierFlows.AbstractDiagnostic}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.increment!",
    "category": "Method",
    "text": "increment!(diag)\n\nIncrement the Diagnostic diag.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.savediagnostic-Tuple{FourierFlows.AbstractDiagnostic,String,String}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.savediagnostic",
    "category": "Method",
    "text": "savediagnostic(diag, diagname)\n\nSave diagnostics to file, labeled by the string diagname.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveoutput-Tuple{FourierFlows.Output}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveoutput",
    "category": "Method",
    "text": "saveoutput(out)\n\nSave current output fields for file in out.filename.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveproblem-Tuple{FourierFlows.AbstractProblem,String}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveproblem",
    "category": "Method",
    "text": "saveproblem(prob, filename)\n\nSave certain aspects of a problem timestepper, grid, and params. Functions that are fields in params are not saved.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.stepforward!-Tuple{FourierFlows.Problem,AbstractArray,Any}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.stepforward!",
    "category": "Method",
    "text": "stepforward!(prob, diags, nsteps)\n\nStep forward the problem prob for nsteps while calculating the  diagnostics in diags.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.stepforward!-Tuple{FourierFlows.Problem,Any}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.stepforward!",
    "category": "Method",
    "text": "stepforward!(prob, nsteps)\n\nStep forward the problem 'prob' for 'nsteps'.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.stepforward!-Tuple{FourierFlows.Problem}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.stepforward!",
    "category": "Method",
    "text": "stepforward!(prob)\n\nStep forward the Problem prob for one timestep.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.update!-Tuple{FourierFlows.AbstractDiagnostic}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.update!",
    "category": "Method",
    "text": "update!(diag)\n\nUpdate diag with its current value.\n\n\n\n"
},

{
    "location": "man/docstrings.html#Functions-exported-from-FourierFlows:-1",
    "page": "Functions exported from FourierFlows:",
    "title": "Functions exported from FourierFlows:",
    "category": "section",
    "text": "Modules = [FourierFlows]\nPrivate = false\nOrder = [:function]"
},

{
    "location": "man/docstrings.html#FourierFlows.ZeroDGrid",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.ZeroDGrid",
    "category": "Type",
    "text": "Doc.\n\n\n\n"
},

{
    "location": "man/docstrings.html#Private-types-in-module-FourierFlows:-1",
    "page": "Functions exported from FourierFlows:",
    "title": "Private types in module FourierFlows:",
    "category": "section",
    "text": "Modules = [FourierFlows]\nPublic = false\nOrder = [:type]"
},

]}
