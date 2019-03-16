module Seek

import Statistics

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

function normalize(a)
    m=[Statistics.median(skipmissing(r)) for r in eachrow(a)]
    abs.(hcat((c-m for c in eachcol(a))...))
end



end # module
