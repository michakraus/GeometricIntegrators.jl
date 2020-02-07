
using Zarr
using Zarr: writechunk!, NoCompressor

const nd = 10
const nt = 1000
const ni = 10

pathname = "test_zarr"

xs = rand(nd,nt,ni)

ds = Zarr.DirectoryStore(pathname)
dx = zcreate(Float64, ds, nd, nt, ni; chunks=(nd,1,1), compressor=NoCompressor())


function copy_to_zarr_1(ds, data)
    ds[:,:,:] = data
end

function copy_to_zarr_2(ds, data)
    for i in axes(data,2)
        ds[:,i,:] = data[:,i,:]
    end
end

function copy_to_zarr_3(ds, data)
    for j in axes(data,3)
        ds[:,:,j] = data[:,:,j]
    end
end

function copy_to_zarr_4(ds, data)
    for i in axes(data,2)
        for j in axes(data,3)
            ds[:,i,j] = data[:,i,j]
        end
    end
end

function copy_to_zarr_5(ds, data)
    tx = zeros(eltype(data), size(data,1))

    for j in axes(data,3)
        for i in axes(data,2)
            for k in axes(data,1)
                tx[k] = data[k,i,j]
            end

            writechunk!(tx, ds, CartesianIndex(1,i,j))
        end
    end
end


#append!(dx,xs,dims=2)
#@time append!(dx,xs,dims=2)

copy_to_zarr_1(dx, xs)
@time copy_to_zarr_1(dx, xs)

copy_to_zarr_2(dx, xs)
@time copy_to_zarr_2(dx, xs)

copy_to_zarr_3(dx, xs)
@time copy_to_zarr_3(dx, xs)

copy_to_zarr_4(dx, xs)
@time copy_to_zarr_4(dx, xs)

copy_to_zarr_5(dx, xs)
@time copy_to_zarr_5(dx, xs)


rm(pathname, recursive=true)



# release 06.02.2020
#   3.280983 seconds (730.02 k allocations: 38.572 MiB, 0.89% gc time)
#   3.231358 seconds (759.49 k allocations: 40.532 MiB, 0.79% gc time)
#   3.671065 seconds (730.30 k allocations: 39.346 MiB, 0.82% gc time)
#   3.588933 seconds (1.02 M allocations: 52.531 MiB, 0.66% gc time)
#   2.658375 seconds (340.00 k allocations: 15.395 MiB)
