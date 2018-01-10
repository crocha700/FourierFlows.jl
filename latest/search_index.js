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
    "text": "FourierFlows is a registered package and so can be installed via Pkg.add.Pkg.add(\"FourierFlows\")For now, this package supports Julia 0.6."
},

{
    "location": "index.html#Source-code-organization-1",
    "page": "Home",
    "title": "Source code organization",
    "category": "section",
    "text": "The code is divided along conceptual lines into problem-agnostic and problem-specific components. Files that contain problem-agnostic parts of the code are stored in /src. Files in /src define the domain, 'AbstractTypes' that supertype problem-specific types, and time-stepper types and routines. Problem-specific modules are stores in /src/physics.Here's an overview of the code structure:/src/\nFourierFlows.jl\nDefines supertyping AbstractParams, AbstractGrid, etc.\nDefines a Problem type to organize the grid, vars, params,   equation, and timestepper into a single structure.\nIncludes all sources files and physics files.\ntimesteppers.jl: defines modules and stepforward! routines for   various time-steppers. Current implemented time-steppers are:\nForward Euler (+ Filtered Forward Euler)\n3rd-order Adams-Bashforth (AB3)\n4th-order Runge-Kutta (RK4)\n4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)\n(+ Filtered ETDRK4)\nphysics/\ntwodturb.jl: Defines a TwoDTurb module that provides a       solver for the two-dimensional vorticity equation.\nbarotropicqg.jl: Defines a BarotropicQG module that provides       several solvers for the barotropic QG model that permit beta,       topography, beta + topography, and forcing.\ntwomodeboussinesq.jl: Defines a TwoModeBoussinesq module       that provides solvers for a two-mode truncation of the       rotating, stratified Boussinesq equation.\nniwqg.jl: Defines a NIWQG module that provides a solver       for the vertical-plane-wave model for the interaction of       a near-inertial wave field and quasi-geostrophic flow.\ntraceradvdiff.jl: Defines a TracerAdvDiff module that       provides a solver for a two-dimensional and periodic tracer       field in a given 2D flow (u, w), which can be an arbitrary       function of x, z, and t."
},

{
    "location": "index.html#Writing-fast-solvers-1",
    "page": "Home",
    "title": "Writing fast solvers",
    "category": "section",
    "text": "The performance-intensive part of the code involves just two functions: the timestepping scheme stepforward!, and the function calcNL! that calculates the nonlinear part of the given equation's right-hand side. Optimization of these two functions for a given problem will produce the fastest possible code."
},

{
    "location": "index.html#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "An example script that simulates decaying two-dimensional turbulence reproducing the results of the paper byMcWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. _J. Fluid Mech._, 146, 21-43.is found in examples/twodturb/McWilliams.jl."
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
    "text": "Resize the diagnostic data and time arrays. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.groupsize-Tuple{JLD2.Group}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.groupsize",
    "category": "Method",
    "text": "Find the number of elements in a JLD2 group. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.increment!-Tuple{AbstractArray}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.increment!",
    "category": "Method",
    "text": "Increment an array of diagnostics. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.increment!-Tuple{FourierFlows.AbstractDiagnostic}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.increment!",
    "category": "Method",
    "text": "Increment a diagnostic by calculating a new value. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.savediagnostic-Tuple{FourierFlows.AbstractDiagnostic,String,String}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.savediagnostic",
    "category": "Method",
    "text": "Save diagnostics to file.\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveoutput-Tuple{AbstractArray}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveoutput",
    "category": "Method",
    "text": "Save an array of outputs to file. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveoutput-Tuple{FourierFlows.OldOutput}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveoutput",
    "category": "Method",
    "text": "Save output to file. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveoutput-Tuple{FourierFlows.Output}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveoutput",
    "category": "Method",
    "text": "Save the current output fields. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveproblem-Tuple{FourierFlows.AbstractProblem,String}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveproblem",
    "category": "Method",
    "text": "Save certain aspects of a Problem. Entire problems cannot be saved in general, because functions cannot be saved (and functions may use arbitrary numbers of global variables that cannot be included in a saved object).\n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.saveproblem-Tuple{FourierFlows.Output}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.saveproblem",
    "category": "Method",
    "text": "Save attributes of the Problem associated with the given Output. \n\n\n\n"
},

{
    "location": "man/docstrings.html#FourierFlows.update!-Tuple{FourierFlows.AbstractDiagnostic}",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.update!",
    "category": "Method",
    "text": "Update the current value of the diagnostic. \n\n\n\n"
},

{
    "location": "man/docstrings.html#Functions-exported-from-FourierFlows:-1",
    "page": "Functions exported from FourierFlows:",
    "title": "Functions exported from FourierFlows:",
    "category": "section",
    "text": "Modules = [FourierFlows]\nPrivate = false\nOrder = [:function]"
},

{
    "location": "man/docstrings.html#FourierFlows.OldOutput",
    "page": "Functions exported from FourierFlows:",
    "title": "FourierFlows.OldOutput",
    "category": "Type",
    "text": "Original output type for FourierFlows problems. \n\n\n\n"
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
