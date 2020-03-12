
using NCDatasets

const nd = 10
const nt = 10000
const ni = 10

file = "test.nc"

xs = rand(nd,nt,ni)

ds = Dataset(file, "c")
dd = defDim(ds,"d",nd)
dt = defDim(ds,"t",nt)
di = defDim(ds,"i",ni)
dx = defVar(ds,"x",Float64,("d","t","i"))


function copy_to_nc_1(ds, data)
    ds[:,:,:] = data
end

function copy_to_nc_2(ds, data)
    for j in axes(data,3)
        ds[:,:,j] = data[:,:,j]
    end
end

function copy_to_nc_3(ds, data)
    for i in axes(data,2)
        ds[:,i,:] = data[:,i,:]
    end
end

function copy_to_nc_4(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            ds[:,i,j] = data[:,i,j]
        end
    end
end

function copy_to_nc_5(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            for k in axes(data,1)
                ds[k,i,j] = data[k,i,j]
            end
        end
    end
end


copy_to_nc_1(dx,xs)
@time copy_to_nc_1(dx,xs)

copy_to_nc_2(dx,xs)
@time copy_to_nc_2(dx,xs)

copy_to_nc_3(dx,xs)
@time copy_to_nc_3(dx,xs)

copy_to_nc_4(dx,xs)
@time copy_to_nc_4(dx,xs)

copy_to_nc_5(dx,xs)
@time copy_to_nc_5(dx,xs)


close(ds)

rm(file)



# release 06.02.2020
#   0.006205 seconds (22 allocations: 7.631 MiB)
#   0.012176 seconds (694 allocations: 15.301 MiB)
#   2.062490 seconds (697.96 k allocations: 58.105 MiB, 0.22% gc time)
#   3.588344 seconds (7.18 M allocations: 455.926 MiB, 1.21% gc time)
#  28.099260 seconds (17.90 M allocations: 975.003 MiB, 0.32% gc time)


# #master 06.02.2020
#   0.015616 seconds (18 allocations: 15.260 MiB, 51.76% gc time)
#   0.012678 seconds (244 allocations: 22.906 MiB)
#   1.691456 seconds (210.00 k allocations: 40.894 MiB, 3.70% gc time)
#   3.096364 seconds (2.40 M allocations: 207.520 MiB, 0.54% gc time)
#  26.199991 seconds (7.95 M allocations: 411.208 MiB, 0.17% gc time)


# #return-type 06.02.2020
#   0.006005 seconds (22 allocations: 7.631 MiB)
#   0.009848 seconds (284 allocations: 15.281 MiB)
#   1.684312 seconds (269.49 k allocations: 37.529 MiB, 0.34% gc time)
#   3.121493 seconds (2.99 M allocations: 244.063 MiB, 0.92% gc time)
#  28.465680 seconds (17.90 M allocations: 975.003 MiB, 0.41% gc time)
