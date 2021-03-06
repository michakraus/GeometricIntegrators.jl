
"createHDF5: Creates or opens HDF5 file."
function createHDF5(sol::Solution, file::AbstractString; overwrite=true)
    if overwrite
        flag = "w"
        get_config(:verbosity) > 1 ? @info("Creating HDF5 file ", file) : nothing
        isfile(file) ? @warn("Overwriting existing HDF5 file.") : nothing
    else
        flag = "r+"
        get_config(:verbosity) > 1 ? @info("Opening HDF5 file ", file) : nothing
    end

    # create or open HDF5 file
    h5 = h5open(file, flag)

    return h5
end

"save_attributes: Saves common attributes of Solution to HDF5 file."
function save_attributes(sol::Solution)
    save_attributes(sol, hdf5(sol))
end

"save_attributes: Saves attributes of Deterministic Solutions to HDF5 file."
function save_attributes(sol::DeterministicSolution, h5::HDF5.File)
    # save attributes
    attributes(h5)["ntime"] = ntime(sol)
    attributes(h5)["nsave"] = nsave(sol)
    attributes(h5)["nsamples"]  = nsamples(sol)
end

"write_to_hdf5: Wrapper for saving Solution to HDF5 file."
function Common.write_to_hdf5(solution::Solution)
    write_to_hdf5(solution, hdf5(solution), offset(solution))
end

"Creates HDF5 file, writes solution to file, and closes file."
function Common.write_to_hdf5(solution::Solution, file::AbstractString)
    h5 = create_hdf5(solution, file)
    write_to_hdf5(solution, h5)
    close(h5)
end


"Creates HDF5 file and initialises datasets for solution object."
function create_hdf5!(solution::Solution, file::AbstractString; kwargs...)
    solution.h5 = create_hdf5(solution, file; kwargs...)
end

"Creates HDF5 file and initialises datasets for deterministic solution object."
function create_hdf5(solution::DeterministicSolution, file::AbstractString)
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, h5)

    # create dataset
    init_timeteps_in_hdf5(solution, h5)
    init_solution_in_hdf5(solution, h5)

    return h5
end

function init_timeteps_in_hdf5(solution::Solution{DT,TT}, h5::HDF5.File) where {DT,TT}
    t = create_dataset(h5, "t", TT, ((solution.nt+1,), (-1,)), chunk=(1,))
    t[1] = timesteps(solution)[0]
end


function init_solution_in_hdf5(solution::SolutionODE{DT,TT,1}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    q[1] = solution.q[0]
end

function init_solution_in_hdf5(solution::SolutionODE{AT,TT,1}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin])
    elaxes = axes(solution.q[begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    q[elaxes..., 1] = solution.q[0]
end

function init_solution_in_hdf5(solution::SolutionODE{DT,TT,2}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    q[1,:] = solution.q[0,:]
end

function init_solution_in_hdf5(solution::SolutionODE{AT,TT,2}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin,begin])
    elaxes = axes(solution.q[begin,begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.q,2)
        q[elaxes..., 1, k] = solution.q[0, k]
    end
end


function init_solution_in_hdf5(solution::SolutionPODE{DT,TT,1}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    p = create_dataset(h5, "p", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    q[1] = solution.q[0]
    p[1] = solution.p[0]
end

function init_solution_in_hdf5(solution::SolutionPODE{AT,TT,1}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin])
    elaxes = axes(solution.q[begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    q[elaxes..., 1] = solution.q[0]

    elsize = size(solution.p[begin])
    elaxes = axes(solution.p[begin])
    p = create_dataset(h5, "p", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    p[elaxes..., 1] = solution.p[0]
end

function init_solution_in_hdf5(solution::SolutionPODE{DT,TT,2}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    p = create_dataset(h5, "p", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    q[1,:] = solution.q[0,:]
    p[1,:] = solution.p[0,:]
end

function init_solution_in_hdf5(solution::SolutionPODE{AT,TT,2}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin,begin])
    elaxes = axes(solution.q[begin,begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.q,2)
        q[elaxes..., 1, k] = solution.q[0, k]
    end

    elsize = size(solution.p[begin,begin])
    elaxes = axes(solution.p[begin,begin])
    p = create_dataset(h5, "p", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.p,2)
        p[elaxes..., 1, k] = solution.p[0, k]
    end
end


function init_solution_in_hdf5(solution::SolutionDAE{DT,TT,1}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    λ = create_dataset(h5, "λ", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    q[1] = solution.q[0]
    λ[1] = solution.λ[0]
end

function init_solution_in_hdf5(solution::SolutionDAE{AT,TT,1}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin])
    elaxes = axes(solution.q[begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    q[elaxes..., 1] = solution.q[0]

    elsize = size(solution.λ[begin])
    elaxes = axes(solution.λ[begin])
    λ = create_dataset(h5, "λ", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    λ[elaxes..., 1] = solution.λ[0]
end

function init_solution_in_hdf5(solution::SolutionDAE{DT,TT,2}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    q[1,:] = solution.q[0,:]
    λ[1,:] = solution.λ[0,:]
end

function init_solution_in_hdf5(solution::SolutionDAE{AT,TT,2}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin,begin])
    elaxes = axes(solution.q[begin,begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.q,2)
        q[elaxes..., 1, k] = solution.q[0, k]
    end

    elsize = size(solution.λ[begin,begin])
    elaxes = axes(solution.λ[begin,begin])
    λ = create_dataset(h5, "λ", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.λ,2)
        λ[elaxes..., 1, k] = solution.λ[0, k]
    end
end


function init_solution_in_hdf5(solution::SolutionPDAE{DT,TT,1}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    p = create_dataset(h5, "p", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    λ = create_dataset(h5, "λ", DT, ((solution.nt+1,), (-1,)), chunk=(1,))
    q[1] = solution.q[0]
    p[1] = solution.p[0]
    λ[1] = solution.λ[0]
end

function init_solution_in_hdf5(solution::SolutionPDAE{AT,TT,1}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin])
    elaxes = axes(solution.q[begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    q[elaxes..., 1] = solution.q[0]

    elsize = size(solution.p[begin])
    elaxes = axes(solution.p[begin])
    p = create_dataset(h5, "p", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    p[elaxes..., 1] = solution.p[0]

    elsize = size(solution.λ[begin])
    elaxes = axes(solution.λ[begin])
    λ = create_dataset(h5, "λ", DT, ((elsize..., solution.nt+1), (elsize..., -1)), chunk=(elsize..., 1))
    λ[elaxes..., 1] = solution.λ[0]
end

function init_solution_in_hdf5(solution::SolutionPDAE{DT,TT,2}, h5::HDF5.File) where {DT <: Number, TT}
    q = create_dataset(h5, "q", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    p = create_dataset(h5, "p", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    λ = create_dataset(h5, "λ", DT, ((solution.nt+1, solution.ni), (-1, solution.ni)), chunk=(1,1))
    q[1,:] = solution.q[0,:]
    p[1,:] = solution.p[0,:]
    λ[1,:] = solution.λ[0,:]
end

function init_solution_in_hdf5(solution::SolutionPDAE{AT,TT,2}, h5::HDF5.File) where {DT, AT <: Array{DT}, TT}
    elsize = size(solution.q[begin,begin])
    elaxes = axes(solution.q[begin,begin])
    q = create_dataset(h5, "q", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.q,2)
        q[elaxes..., 1, k] = solution.q[0, k]
    end

    elsize = size(solution.p[begin,begin])
    elaxes = axes(solution.p[begin,begin])
    p = create_dataset(h5, "p", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.p,2)
        p[elaxes..., 1, k] = solution.p[0, k]
    end

    elsize = size(solution.λ[begin,begin])
    elaxes = axes(solution.λ[begin,begin])
    λ = create_dataset(h5, "λ", DT, ((elsize..., solution.nt+1, solution.ni), (elsize..., -1, solution.ni)), chunk=(elsize..., 1, 1))
    for k in axes(solution.λ,2)
        λ[elaxes..., 1, k] = solution.λ[0, k]
    end
end


function copy_timeteps_to_hdf5(sol::Solution, h5::HDF5.File, j1, j2, n1, n2)
    if size(h5["t"],1) < j2
        HDF5.set_extent_dims(h5["t"], (j2,))
    end
    h5["t"][j1:j2] = timesteps(sol)[n1:n2]
end


function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,1},SolutionDAE{DT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (j2,))
    end
    h5["q"][j1:j2] = solution.q[n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{AT,TT,1},SolutionDAE{AT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin])
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["q"][elaxes..., j] = solution.q[n]
    end
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{DT,TT,2},SolutionDAE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"],1) < j2
        HDF5.set_extent_dims(h5["q"], (j2, size(h5["q"],2)))
    end
    h5["q"][j1:j2, :] = solution.q[n1:n2, :]
end

function copy_solution_to_hdf5(solution::Union{SolutionODE{AT,TT,2},SolutionDAE{AT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin,begin])
    if size(h5["q"])[end-1] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[begin:end-2]..., j2, size(h5["q"])[end]))
    end
    for k in 1:nsamples(solution.q)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["q"][elaxes..., j, k] = solution.q[n, k]
        end
    end
end


function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,1},SolutionPDAE{DT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (j2,))
    end
    if size(h5["p"])[end] < j2
        HDF5.set_extent_dims(h5["p"], (j2,))
    end
    h5["q"][j1:j2] = solution.q[n1:n2]
    h5["p"][j1:j2] = solution.p[n1:n2]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{AT,TT,1},SolutionPDAE{AT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin])
    if size(h5["q"])[end] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["q"][elaxes..., j] = solution.q[n]
    end

    elaxes = axes(solution.p[begin])
    if size(h5["p"])[end] < j2
        HDF5.set_extent_dims(h5["p"], (size(h5["p"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["p"][elaxes..., j] = solution.p[n]
    end
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{DT,TT,2},SolutionPDAE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["p"],1) < j2
        HDF5.set_extent_dims(h5["p"], (j2, size(h5["p"],2)))
    end
    if size(h5["p"],1) < j2
        HDF5.set_extent_dims(h5["p"], (j2, size(h5["p"],2)))
    end
    h5["q"][j1:j2, :] = solution.q[n1:n2, :]
    h5["p"][j1:j2, :] = solution.p[n1:n2, :]
end

function copy_solution_to_hdf5(solution::Union{SolutionPODE{AT,TT,2},SolutionPDAE{AT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.q[begin,begin])
    if size(h5["q"])[end-1] < j2
        HDF5.set_extent_dims(h5["q"], (size(h5["q"])[begin:end-2]..., j2, size(h5["q"])[end]))
    end
    for k in 1:nsamples(solution.q)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["q"][elaxes..., j, k] = solution.q[n, k]
        end
    end

    elaxes = axes(solution.p[begin,begin])
    if size(h5["p"])[end-1] < j2
        HDF5.set_extent_dims(h5["p"], (size(h5["p"])[begin:end-2]..., j2, size(h5["p"])[end]))
    end
    for k in 1:nsamples(solution.p)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["p"][elaxes..., j, k] = solution.p[n, k]
        end
    end
end


function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,1},SolutionPDAE{DT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["λ"])[end] < j2
        HDF5.set_extent_dims(h5["λ"], (j2,))
    end
    h5["λ"][j1:j2] = solution.λ[n1:n2]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{AT,TT,1},SolutionPDAE{AT,TT,1}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.λ[begin])
    if size(h5["λ"])[end] < j2
        HDF5.set_extent_dims(h5["λ"], (size(h5["λ"])[1:end-1]..., j2))
    end
    for i in eachindex(j1:j2, n1:n2)
        j = (j1:j2)[i]
        n = (n1:n2)[i]
        h5["λ"][elaxes..., j] = solution.λ[n]
    end
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{DT,TT,2},SolutionPDAE{DT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT <: Number, TT}
    if size(h5["λ"],1) < j2
        HDF5.set_extent_dims(h5["λ"], (j2, size(h5["λ"],2)))
    end
    h5["λ"][j1:j2, :] = solution.λ[n1:n2, :]
end

function copy_multiplier_to_hdf5(solution::Union{SolutionDAE{AT,TT,2},SolutionPDAE{AT,TT,2}}, h5::HDF5.File, j1, j2, n1, n2) where {DT, AT <: Array{DT}, TT}
    elaxes = axes(solution.λ[begin,begin])
    if size(h5["λ"])[end-1] < j2
        HDF5.set_extent_dims(h5["λ"], (size(h5["λ"])[begin:end-2]..., j2, size(h5["λ"])[end]))
    end
    for k in 1:nsamples(solution.λ)
        for i in eachindex(j1:j2, n1:n2)
            j = (j1:j2)[i]
            n = (n1:n2)[i]
            h5["λ"][elaxes..., j, k] = solution.λ[n, k]
        end
    end
end


"Append solution to HDF5 file."
function Common.write_to_hdf5(solution::DeterministicSolution, h5::HDF5.File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end

"Append solution to HDF5 file."
function Common.write_to_hdf5(solution::Union{SolutionDAE,SolutionPDAE}, h5::HDF5.File=hdf5(solution), offset=offset(solution))
    # set convenience variables and compute ranges
    j1 = offset+2
    j2 = offset+1+solution.nt

    # copy data from solution to HDF5 dataset
    copy_timeteps_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_solution_to_hdf5(solution, h5, j1, j2, 1, solution.nt)
    copy_multiplier_to_hdf5(solution, h5, j1, j2, 1, solution.nt)

    return nothing
end
