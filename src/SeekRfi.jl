module SeekRfi

import Statistics
import Images

#import FITSIO

#function write_fits(T::DataType, v, name::String)
#    f=FITSIO.FITS(name, "w")
#    write(f, convert(Array{T},transpose(v)))
#    close(f)
#end


function gaussian_filter(V, mask; kernel_m=40, kernel_n=20, sigma_m=7.5, sigma_n=15)
    wd(n,m;sigma_n, sigma_m)=exp.(-n.^2/(2*sigma_n^2).-m.^2/(2*sigma_m^2))

    Vp=zeros(eltype(V), size(V,1)+kernel_n, size(V,2)+kernel_m)
    Vp[kernel_n÷2+1:kernel_n÷2+size(V,1),kernel_m÷2+1:kernel_m÷2+size(V,2)] = V[:,:]
    
    Wfp=zeros(eltype(V), size(V,1)+kernel_n, size(V,2)+kernel_m)
    Wfp[kernel_n÷2+1:kernel_n÷2+size(V,1),kernel_m÷2+1:kernel_m÷2+size(V,2)] = (x->if x 1.0 else 0.0 end).((~).(mask))
    
    Vh=zeros(eltype(V), size(V,1)+kernel_n, size(V,2)+kernel_m)
    Vh2=zeros(eltype(V), size(V,1)+kernel_n, size(V,2)+kernel_m)
    n=collect(-kernel_n÷2:kernel_n÷2)
    m=collect(-kernel_m÷2:kernel_m÷2)
    #println(size(n))

    kernel_0=wd(n, 0, sigma_n=sigma_n, sigma_m=sigma_m)
    kernel_1=wd(0, m, sigma_n=sigma_n, sigma_m=sigma_m)
    
    Vh=_gaussian_filter(Vp, size(V, 1), size(V, 2), Wfp, mask, Vh, Vh2, kernel_0, kernel_1, kernel_m, kernel_n)
    
    Vh=Vh[kernel_n÷2+1:kernel_n÷2+size(V,1),kernel_m÷2+1:kernel_m÷2+size(V,2)]
    Vh[mask]=V[mask]
    Vh    
end


function _gaussian_filter(Vp, vs0, vs1, Wfp, mask, Vh, Vh2, kernel_0, kernel_1, M, N)
    n2=N÷2
    m2=M÷2
    for i in N÷2+1:vs0+N÷2
        for j in M÷2+1:vs1+M÷2
            if mask[i-n2, j-m2]
                Vh[i,j]=zero(eltype(Vh))
            else
                val=sum(Wfp[i-n2:i+n2, j].*Vp[i-n2:i+n2, j].*kernel_0)
                #print(sum(Vp[i-n2+1:i+n2, j])," ")
                
                Vh[i,j]=val/sum(Wfp[i-n2:i+n2, j].*kernel_0)
                #print(Vh[i,j]," ")
            end
        end
    end
    for j2 in M÷2+1:vs1+M÷2
        for i2 in N÷2+1:vs0+N÷2
            if mask[i2-n2, j2-m2]
                Vh2[i2,j2]=zero(eltype(Vh))
            else
                val=sum(Wfp[i2, j2-m2:j2+m2].*Vh[i2, j2-m2:j2+m2].*kernel_1)
                Vh2[i2,j2]=val/sum(Wfp[i2, j2-m2:j2+m2].*kernel_1)
            end
        end
    end
    Vh2
end

function create_masked_array(a::Array{T}, m)::Array{Union{Missing, T}} where {T}
    b=convert(Array{Union{Missing, T}}, a)
    b[m].=missing
    b
end

function strip_mask(a::Array{Union{Missing, T}}; default=zero(T))::Tuple{Array{T}, BitArray} where {T}
    m=ismissing.(a)
    a[m].=default
    (convert(Array{T}, a), m)
end

function normalize(a, mask)
    ma=create_masked_array(a, mask)
    med=[Statistics.median(skipmissing(r)) for r in eachrow(ma)]
    c=abs.(hcat((c-med for c in eachcol(ma))...))
    c=convert(Array{Union{Missing, eltype(a)}}, c)
    #println(typeof(c))
    c,_=strip_mask(c)
    c[mask]=a[mask]
    c
end

function binary_mask_dilation(mask, kernel_m::Int, kernel_n::Int)::BitArray{2}
    #d=Images.dilate(mask)-mask
    pt=ones(kernel_m, kernel_n)
    #(mask+Images.imfilter(d,Images.centered(pt))).>0
    result=Images.imfilter(mask, Images.centered(pt)).>0
    result=Images.imfilter(result, Images.centered(pt)).>0
end

function _sumthreshold(data, mask::BitArray{2}, i::Integer, chi)::BitArray{2}
    tmp_mask=copy(mask)
    ds0, ds1=size(data)
    T=eltype(data)
    for x in 0:ds0-1
        sum = zero(T)
        cnt = 0
        
        for ii in 0:i-1
            if mask[x+1, ii+1] != true
                sum += data[x+1, ii+1]
                cnt += 1
            end
        end
        
        for y in i:ds1-1
            if sum > chi * cnt
                for ii2 in 0:i-1
                    tmp_mask[x+1, y-ii2] = true
                end
            end
                    
            if mask[x+1, y+1] != true
                sum += data[x+1, y+1]
                cnt += 1
            end
            
            if mask[x+1, y-i+1] != true
                sum -= data[x+1, y-i+1]
                cnt -= 1
            end
        end
    end
    tmp_mask
end

function _run_sumthreshold(data, init_mask, eta, M, chi_i; kernel_m=40,kernel_n=20, sigma_m=7.5, sigma_n=15.0)
    smoothed_data=gaussian_filter(data, init_mask;kernel_m=kernel_m, kernel_n=kernel_n, sigma_m=sigma_m, sigma_n=sigma_n)
    res=data-smoothed_data

    st_mask=copy(init_mask)

    for (m, chi) in zip(M, chi_i)
        #println("m= ",m)
        chi = chi/eta
        if m==1
            st_mask = st_mask .| (res .>= chi)
            #write_fits(Int, st_mask, "mask_0.fits")
            #exit()
        else
            #println(size(res))
            #println(size(st_mask))
            st_mask = _sumthreshold(res, st_mask, m, chi)
            st_mask = convert(BitArray, 
            transpose(_sumthreshold(convert(Array, transpose(res)), convert(BitArray, transpose(st_mask)), m, 
            chi))
            )
            write_fits(Int, st_mask, "mask_"*string(m)*".fits")
        end
    end
    write_fits(Int, st_mask, "final_"*string(eta)*".fits")
    st_mask
end

function get_rfi_mask(data; mask::Union{Missing, BitArray{2}}=missing,kernel_m=40, kernel_n=20, sigma_m=7.5, sigma_n=15, di_args=(3,7), chi_1=35000.0, eta_i=[0.5, 0.55, 0.62, 0.75, 1], 
    normalize_standing_waves=true, suppress_dilation=false)
    T=eltype(data)

    if ismissing(mask)
        mask=falses(size(data)...)
    end

    if normalize_standing_waves
        data=normalize(data, mask)
    end
    p=convert(T,1.5)
    m=1:7
    M=2 .^ (m .-1)
    chi_i = chi_1 ./ p .^ log2.(m)
    st_mask=copy(mask)
    for (i,eta) in enumerate(eta_i)
        st_mask=_run_sumthreshold(data, st_mask, eta, M , chi_i; kernel_m=kernel_m,kernel_n=kernel_n, sigma_m=sigma_m, sigma_n=sigma_n)
        #write_fits(Int, st_mask, "final_"*string(eta)*".fits")
    end

    if suppress_dilation
        st_mask
    else
        dilated_mask=st_mask
        dilated_mask = binary_mask_dilation(dilated_mask .⊻ mask, di_args...)

        dilated_mask .| st_mask
    end
end


end # module
