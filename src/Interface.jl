using LinearAlgebra
using DiffEqBase
using LsqFit

abstract type AbstractInterface end
export AbstractInterface

abstract type AbstractBoundaryInterface <: AbstractInterface end
export AbstractBoundaryInterface

abstract type AbstractInternalInterface <: AbstractInterface end
export AbstractInternalInterface

abstract type AbstractReactiveInternalInterface <: AbstractInternalInterface end
export AbstractReactiveInternalInterface

struct EmptyInterface <: AbstractInterface end
export EmptyInterface

struct ReactiveInternalInterface{T,B,C,C2,N,Q<:AbstractReaction,X} <: AbstractReactiveInternalInterface
    domain1::T
    domain2::N
    reactions::Array{Q,1}
    veckinetics::X
    veckineticsinds::Array{Int64,1}
    rxnarray::B
    stoichmatrix::C
    Nrp1::C2
    Nrp2::C2
    A::Float64
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
    reversibililty::Array{Bool,1}
end
function ReactiveInternalInterface(domain1,domain2,reactions,A)
    vectuple,vecinds,otherrxns,otherrxninds,posinds = getveckinetics(reactions)
    rxns = vcat(reactions[vecinds],reactions[otherrxninds])
    rxns = [ElementaryReaction(index=i,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
        productinds=rxn.productinds,kinetics=rxn.kinetics,radicalchange=rxn.radicalchange,reversible=rxn.reversible,pairs=rxn.pairs) for (i,rxn) in enumerate(rxns)]
    rxnarray = getinterfacereactioninds(domain1,domain2,rxns)
    M,Nrp1,Nrp2 = getstoichmatrix(domain1,domain2,reactions)
    reversibility = Array{Bool,1}(getfield.(rxns,:reversible))
    return ReactiveInternalInterface(domain1,domain2,
            rxns,vectuple,posinds,rxnarray,M,Nrp1,Nrp2,A,[1,length(reactions)],
            [0,1],ones(length(rxns)),reversibility),ones(length(rxns))
end
export ReactiveInternalInterface

function getkfskrevs(ri::ReactiveInternalInterface,T1,T2,phi1,phi2,Gs1,Gs2,cstot::Array{Q,1}) where {Q}
    kfs = getkfs(ri,T1,0.0,0.0,Array{Q,1}(),ri.A,phi1)
    Kc = getKc.(ri.reactions,ri.domain1.phase,ri.domain2.phase,Ref(Gs1),Ref(Gs2),T1,phi1)
    krevs = kfs./Kc
    return kfs,krevs
end

function evaluate(ri::ReactiveInternalInterface,dydt,domains,T1,T2,phi1,phi2,Gs1,Gs2,cstot,p::W) where {W<:DiffEqBase.NullParameters}
    kfs,krevs = getkfskrevs(ri,T1,T2,phi1,phi2,Gs1,Gs2,cstot)
    addreactionratecontributions!(dydt,ri.rxnarray,cstot,kfs,krevs,ri.A)
end

function evaluate(ri::ReactiveInternalInterface,dydt,domains,T1,T2,phi1,phi2,Gs1,Gs2,cstot,p)
    kfs,krevs = getkfskrevs(ri,T1,T2,phi1,phi2,Gs1,Gs2,cstot)
    addreactionratecontributions!(dydt,ri.rxnarray,cstot,kfs.*p[ri.parameterindexes[1]:ri.parameterindexes[2]],krevs.*p[ri.parameterindexes[1]:ri.parameterindexes[2]],ri.A)
end
export evaluate


struct ReactiveInternalInterfaceConstantTPhi{J,N,B,B2,B3,C,C2,Q<:AbstractReaction} <: AbstractReactiveInternalInterface
    domain1::J
    domain2::N
    reactions::Array{Q,1}
    rxnarray::B
    stoichmatrix::C
    Nrp1::C2
    Nrp2::C2
    kfs::B2
    krevs::B3
    T::Float64
    A::Float64
    parameterindexes::Array{Int64,1}
    domaininds::Array{Int64,1}
    p::Array{Float64,1}
    reversibility::Array{Bool,1}
end
function ReactiveInternalInterfaceConstantTPhi(domain1,domain2,reactions,T,A,phi=0.0)
    @assert domain1.T == domain2.T 
    reactions = upgradekinetics(reactions,domain1,domain2)
    rxnarray = getinterfacereactioninds(domain1,domain2,reactions)
    kfs = getkf.(reactions,nothing,T,0.0,0.0,Ref([]),A,phi)
    Kc = getKc.(reactions,domain1.phase,domain2.phase,Ref(domain1.Gs),Ref(domain2.Gs),T,phi)
    krevs = kfs./Kc
    M,Nrp1,Nrp2 = getstoichmatrix(domain1,domain2,reactions)
    reversibility = Array{Bool,1}(getfield.(reactions,:reversible))
    return ReactiveInternalInterfaceConstantTPhi(domain1,domain2,reactions,
            rxnarray,M,Nrp1,Nrp2,kfs,krevs,T,A,[1,length(reactions)],
            [0,1],kfs[1:end],reversibility),kfs[1:end]
end
export ReactiveInternalInterfaceConstantTPhi

function getkfskrevs(ri::ReactiveInternalInterfaceConstantTPhi,T1,T2,phi1,phi2,Gs1,Gs2,cstot)
    return ri.kfs,ri.krevs
end

function evaluate(ri::ReactiveInternalInterfaceConstantTPhi,dydt,domains,T1,T2,phi1,phi2,Gs1,Gs2,cstot,p::W) where {W<:DiffEqBase.NullParameters}
    addreactionratecontributions!(dydt,ri.rxnarray,cstot,ri.kfs,ri.krevs,ri.A)
end

function evaluate(ri::ReactiveInternalInterfaceConstantTPhi,dydt,domains,T1,T2,phi1,phi2,Gs1,Gs2,cstot,p)
    if all(p[ri.parameterindexes[1]:ri.parameterindexes[2]] .== ri.kfs)
        kfs = ri.kfs
    else 
        kfs = p[ri.parameterindexes[1]:ri.parameterindexes[2]]
    end
    if length(Gs1) == 0 || length(Gs2) == 0 || (all(Gs1 .== ri.domain1.Gs) && all(Gs2 .== ri.domain2.Gs))
        krevs = ri.krevs
    else
        Kc = getKcs(ri,T1,Gs1,Gs2)
        krevs = kfs./Kc
    end
    addreactionratecontributions!(dydt,ri.rxnarray,cstot,kfs,krevs,ri.A)
end
export evaluate

"""
construct the stochiometric matrix for the reactions crossing both domains and the reaction molecule # change
"""
function getstoichmatrix(domain1,domain2,rxns)
    M = spzeros(length(rxns),length(domain1.phase.species)+length(domain2.phase.species))
    Nrp1 = zeros(length(rxns))
    Nrp2 = zeros(length(rxns))
    N1 = length(domain1.phase.species)
    spcs1 = domain1.phase.species
    spcs2 = domain2.phase.species
    for (i,rxn) in enumerate(rxns)
        Nrp1[i] = Float64(length([x for x in rxn.products if x in spcs1]) - length([x for x in rxn.reactants if x in spcs1]))
        Nrp2[i] = Float64(length([x for x in rxn.products if x in spcs2]) - length([x for x in rxn.reactants if x in spcs2]))
        for (j,r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            M[i,isfirst ? ind : ind+N1] += 1
        end
        for (j,r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            M[i,isfirst ? ind : ind+N1] -= 1
        end
    end
    return M,Nrp1,Nrp2
end

function getinterfacereactioninds(domain1,domain2,reactions)
    indices = zeros(Int64,(6,length(reactions)))
    N1 = length(domain1.phase.species)
    for (i,rxn) in enumerate(reactions)
        for (j,r) in enumerate(rxn.reactants)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            indices[j,i] = isfirst ? ind : ind+N1
        end
        for (j,r) in enumerate(rxn.products)
            isfirst = true
            ind = findfirst(isequal(r),domain1.phase.species)
            if ind === nothing
                isfirst = false
                ind = findfirst(isequal(r),domain2.phase.species)
            end
            indices[j+3,i] = isfirst ? ind : ind+N1
        end
    end
    return indices
end

function upgradekinetics(rxns,domain1,domain2)
    domain1surf = hasproperty(domain1.phase,:sitedensity)
    domain2surf = hasproperty(domain2.phase,:sitedensity)
    @assert !(domain1surf && domain2surf)
    if domain1surf
        surfdomain = domain1
    elseif domain2surf
        surfdomain = domain2
    end
    for (i,rxn) in enumerate(rxns)
        if isa(rxn.kinetics,StickingCoefficient)
            spc = [spc for spc in rxn.reactants if !(spc in surfdomain.phase.species)]
            @assert length(spc) == 1
            kin = stickingcoefficient2arrhenius(rxn.kinetics,surfdomain.phase.sitedensity,length(rxn.reactants)-1,spc[1].molecularweight)
            rxns[i] = ElementaryReaction(index=rxn.index,reactants=rxn.reactants,reactantinds=rxn.reactantinds,products=rxn.products,
                productinds=rxn.productinds,kinetics=kin,radicalchange=rxn.radicalchange,reversible=rxn.reversible,pairs=rxn.pairs)
        end
    end
    return rxns
end

function stickingcoefficient2arrhenius(sc,sitedensity,N,mw;Tmin=300.0,Tmax=2000.0)
    mass = mw/Na
    ksc(T) = sc(T)/sitedensity^N*sqrt(kB*T/(2.0*pi*mass))
    Ts = Array{Float64,1}(Tmin:10:Tmax);
    kscvals = ksc.(Ts)
    k(T,p) = log.(abs(p[1]).*T.^p[2].*exp.(-p[3]./(R.*T)))
    p0 = [sc.A/sitedensity*sqrt(kB*1000.0/(2.0*pi*mass)),0.5,sc.Ea]
    fit = curve_fit(k,Ts,log.(kscvals),p0;x_tol=1e-18)
    @assert fit.converged
    p = fit.param
    p[1] = abs(p[1])
    return Arrhenius(;A=p[1],n=p[2],Ea=p[3])
end

struct Inlet{Q<:Real,S,V<:AbstractArray,U<:Real,X<:Real} <: AbstractBoundaryInterface
    domain::S
    y::V
    F::Function
    T::U
    P::X
    H::Q
end

function Inlet(domain::V,conddict::Dict{String,X},F::Function) where {V,X<:Real,B<:Real}
    y = makespcsvector(domain.phase,conddict)
    T = conddict["T"]
    P = conddict["P"]
    yout = y./sum(y)
    H = dot(getEnthalpy.(getfield.(domain.phase.species,:thermo),T),yout)
    return Inlet(domain,yout,F,T,P,H)
end

export Inlet

struct Outlet{V} <: AbstractBoundaryInterface
    domain::V
    F::Function
end
export Outlet


struct Condensation{S,V<:AbstractArray,U<:Real,X<:Real} <: AbstractBoundaryInterface
    domain::S
    x_vap::V
    T_vap::U
    P_vap::X
    kLA::Array{Float64,1}
    kH::Array{Float64,1}
end

function Condensation(domain::V,conddict::Dict{String,X},kLA::Array{Float64,1},kH::Array{Float64,1}) where {V,X<:Real}
    y = makespcsvector(domain.phase,conddict)
    T_vap = conddict["T"]
    P_vap = conddict["P"]
    x_vap = y./sum(y)
    return Condensation(domain,x_vap,T_vap,P_vap,kLA,kH)
end

export Condensation

struct Evaporation{V} <: AbstractBoundaryInterface
    domain::V
    kLA::Array{Float64,1}
    kH::Array{Float64,1}
end

export Evaporation


