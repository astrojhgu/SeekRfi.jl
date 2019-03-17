module SeekRfi

import Statistics
import Images

function gaussian_filter(V, mask; M=40, N=20, sigma_m=0.5, sigma_n=0.5)
    wd(n,m;sigma_n, sigma_m)=exp.(-n.^2/(2*sigma_n^2).-m.^2/(2*sigma_m^2))

    Vp=zeros(eltype(V), size(V,1)+N, size(V,2)+M)
    Vp[N÷2+1:N÷2+size(V,1),M÷2+1:M÷2+size(V,2)] = V[:,:]
    
    Wfp=zeros(eltype(V), size(V,1)+N, size(V,2)+M)
    Wfp[N÷2+1:N÷2+size(V,1),M÷2+1:M÷2+size(V,2)] = (x->if x 1.0 else 0.0 end).((~).(mask))
    
    Vh=zeros(eltype(V), size(V,1)+N, size(V,2)+M)
    Vh2=zeros(eltype(V), size(V,1)+N, size(V,2)+M)
    n=collect(-N÷2+1:N÷2)
    m=collect(-M÷2+1:M÷2)

    kernel_0=wd(n, 0, sigma_n=sigma_n, sigma_m=sigma_m)
    kernel_1=wd(0, m, sigma_n=sigma_n, sigma_m=sigma_m)
    
    Vh=_gaussian_filter(Vp, size(V, 1), size(V, 2), Wfp, mask, Vh, Vh2, kernel_0, kernel_1, M, N)
    
    Vh=Vh[N÷2+1:N÷2+size(V,1),M÷2+1:M÷2+size(V,2)]
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
                val=sum(Wfp[i-n2+1:i+n2, j].*Vp[i-n2+1:i+n2, j].*kernel_0)
                #print(sum(Vp[i-n2+1:i+n2, j])," ")
                
                Vh[i,j]=val/sum(Wfp[i-n2+1:i+n2, j].*kernel_0)
                print(Vh[i,j]," ")
            end
        end
    end
    for j2 in M÷2+1:vs1+M÷2
        for i2 in N÷2+1:vs0+N÷2
            if mask[i2-n2, j2-m2]
                Vh2[i2,j2]=zero(eltype(Vh))
            else
                val=sum(Wfp[i2, j2-m2+1:j2+m2].*Vh[i2, j2-m2+1:j2+m2].*kernel_1)
                Vh2[i2,j2]=val/sum(Wfp[i2, j2-m2+1:j2+m2].*kernel_1)
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

function strip_mask(a::Array{Union{Missing, T}}; default=zero(T))::Tuple{Array{T}, Array{Bool}} where {T}
    m=ismissing.(a)
    a[m].=default
    (convert(Array{T}, a), m)
end

function normalize(a, mask)
    ma=create_masked_array(a, mask)
    med=[Statistics.median(skipmissing(r)) for r in eachrow(ma)]
    c=abs.(hcat((c-med for c in eachcol(ma))...))
    c,_=strip_mask(c)
    c[mask]=a[mask]
    c
end

function binary_mask_dilation(mask::Array{Bool, 2}, ss1::Int, ss2::Int)::Array{Bool, 2}
    #d=Images.dilate(mask)-mask
    pt=trues(ss1, ss2)
    #(mask+Images.imfilter(d,Images.centered(pt))).>0
    Images.imfilter(mask, Images.centered(pt))
end

function _sumthreshold(data, mask::Array{Bool,2}, i, chi, ds0, ds1)::Array{Bool, 2}
    tmp_mask=copy(mask)
    T=eltype(data)
    for x in 1:ds0
        sum=zero(T)
        cnt=0

        for ii in 1:i
            if  ! mask[x, ii]
                sum+=data[x,ii]
                cnt+=1
            end
        end

        for y in i+1:ds1
            if sum > chi*cnt
                for ii2 in 1:i
                    tmp_mask[x, y-ii2]=true
                end
            end

            if ! mask[x,y]
                sum+=data[x,y]
                cnt+=1
            end

            if ! mask[x, y-i]
                sum -= data[x, y-i]
                cnt-=1
            end
        end
    end
    tmp_mask
end

function _run_sumthreshold(data, init_mask, eta, M, chi_i)
    smoothed_data=gaussian_filter(data, init_mask;M=40, N=20, sigma_m=7.5, sigma_n=15)
    res=data-smoothed_data

    st_mask=copy(init_mask)

    for (m, chi) in zip(M, chi_i)
        chi = chi/eta
        if m==1
            st_mask = st_mask .| (res .>= chi)
        else
            st_mask = _sumthreshold(res, st_mask, m, chi, size(res,1), size(res,2))
            st_mask = _sumthreshold(convert(Array, transpose(res)), convert(Array, transpose(st_mask)), 
            chi, size(res, 2), size(res, 1))
        end
    end
    st_mask
end

function get_rfi_mask(data; mask::Union{Missing, Array{Bool, 2}}=missing, ss1=3, ss2=7, chi_1=35000.0, eta_i=[0.5, 0.55, 0.62, 0.75, 1], 
    normalize_standing_waves=true, suppress_dilation=false)
    if ismissing(mask)
        mask=falses(size(data)...)
    end

    if normalize_standing_waves
        data=normalize(data, mask)
    end

    p=convert(T, 1.5)
    m=1:7
    M=convert(T, 2) .^ m
    chi_i = chi_1 / p .^ log2.(m)
    st_mask=copy(mask)
    for eta in eta_i
        st_mask=_run_sumthreshold(data, st_mask, eta, M, chi_i)
    end

    dilated_mask=st_mask
    dilated_mask = binary_mask_dilation(st_mask .⊻ mask, ss1, ss2)

    dilated_mask .| st_mask
end

end # module
