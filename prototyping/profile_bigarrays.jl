
using BigArrays
using BigArrays.BinDicts
using BigArrays.Infos

const nd = 10
const nt = 10000
const ni = 10

pathname = "test_ba"
mkdir(pathname)
mkdir(joinpath(pathname, "x"))

infoString = """
{"num_channels": 1, "type": "image", "data_type": "float64", "scales": [
{"encoding": "gzip", "key": "x", "chunk_sizes": [[10, 1, 1]], "size": [10, 1000, 10], "resolution": [1, 1, 1], "voxel_offset": [0, 0, 0]}
]}
"""

write( joinpath(pathname, "info"), infoString )

ba = BigArray( pathname )


xs = rand(nd,nt,ni)

ba[1:nd,1:nt,1:ni] = xs


function copy_to_ba_1(ds, data)
    ds[:,:,:] = data
end

function copy_to_ba_2(ds, data)
    for j in axes(data,3)
        ds[:,:,j] = data[:,:,j]
    end
end

function copy_to_ba_3(ds, data)
    for i in axes(data,2)
        ds[:,i,:] = data[:,i,:]
    end
end

function copy_to_ba_4(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            ds[:,i,j] = data[:,i,j]
        end
    end
end

function copy_to_ba_5(ds, data)
    for j in axes(data,3)
        for i in axes(data,2)
            for k in axes(data,1)
                ds[k,i,j] = data[k,i,j]
            end
        end
    end
end


# copy_to_ba_1(ba,xs)
# @time copy_to_ba_1(ba,xs)
#
# copy_to_ba_2(ba,xs)
# @time copy_to_ba_2(ba,xs)
#
# copy_to_ba_3(ba,xs)
# @time copy_to_ba_3(ba,xs)
#
# copy_to_ba_4(ba,xs)
# @time copy_to_ba_4(ba,xs)
#
# copy_to_ba_5(ba,xs)
# @time copy_to_ba_5(ba,xs)


# close(ba)

rm(pathname, recursive=true)
