
using HDF5

const nd = 10
const nt = 10000
const ni = 10

file = "test.h5"

xs = rand(nd,nt,ni)

h5 = h5open(file, "w")
dx = d_create(h5, "x", eltype(xs), ((nd,nt,ni),(nd,-1,ni)), "chunk", (nd,1,1))


function copy_to_hdf5_1(ds, data)
    ds[:,:,:] = data
end

function copy_to_hdf5_2(ds, data)
    for j in axes(data,3)
        ds[:,:,j] = data[:,:,j]
    end
end

function copy_to_hdf5_3(ds, data)
    for i in axes(data,2)
        ds[:,i,:] = data[:,i,:]
    end
end

function copy_to_hdf5_4(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            ds[:,i,j] = data[:,i,j]
        end
    end
end

function copy_to_hdf5_5(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            for k in axes(data,1)
                ds[k,i,j] = data[k,i,j]
            end
        end
    end
end


copy_to_hdf5_1(dx,xs)
@time copy_to_hdf5_1(dx,xs)

copy_to_hdf5_2(dx,xs)
@time copy_to_hdf5_2(dx,xs)

copy_to_hdf5_3(dx,xs)
@time copy_to_hdf5_3(dx,xs)

copy_to_hdf5_4(dx,xs)
@time copy_to_hdf5_4(dx,xs)

copy_to_hdf5_5(dx,xs)
@time copy_to_hdf5_5(dx,xs)


close(h5)

rm(file)



# release 07.02.2020
#   1.863731 seconds (287 allocations: 14.000 KiB)
#   1.590241 seconds (2.08 k allocations: 7.728 MiB)
#   1.627206 seconds (2.06 M allocations: 106.338 MiB, 0.57% gc time)
#   5.331719 seconds (13.39 M allocations: 642.239 MiB, 1.62% gc time)
#  32.912910 seconds (84.85 M allocations: 3.976 GiB, 2.48% gc time)
