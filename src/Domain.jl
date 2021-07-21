using Parameters
using LinearAlgebra
using StaticArrays
using Calculus
using DiffEqBase
using ForwardDiff
using Tracker
using ReverseDiff

abstract type AbstractDomain end
export AbstractDomain

abstract type AbstractConstantKDomain <: AbstractDomain  end
export AbstractConstantKDomain

abstract type AbstractVariableKDomain <: AbstractDomain end
export AbstractVariableKDomain

mutable struct ConstantTPDomain{N<:AbstractPhase,S<:Integer,W<:Real, W2<:Real, I<:Integer, Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    P::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end

function ConstantTPDomain(;phase::E2,initialconds::Dict{X,X2},constantspecies::Array{X3,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {E<:Real,E2<:AbstractPhase,Q<:AbstractInterface,W<:Real,X,X2,X3}
    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    y0 = zeros(length(phase.species)+1)
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end


    @assert T != 0.0
    @assert P != 0.0
    ns = y0[1:end-1]
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase,T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = getfield.(phase.species,:diffusion)(T=T,mu=mu,P=P)::Array{typeof(T),1}
    else
        diffs = Array{typeof(T),1}()
    end
    C = P/(R*T)
    V = N*R*T/P
    y0[end] = V
    kfs,krevs = getkfkrevs(phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    kfsp = deepcopy(kfs)
    for ind in efficiencyinds
        kfsp[ind] = 1.0
    end
    p = vcat(deepcopy(Gs),kfsp)
    if sparse
        jacobian=spzeros(typeof(T),length(phase.species),length(phase.species))
    else
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTPDomain(phase,[1,length(phase.species),length(phase.species)+1],[1,length(phase.species)+length(phase.reactions)],constspcinds,
        T,P,kfs,krevs,efficiencyinds,Gs,rxnarray,mu,diffs,jacobian,sensitivity,false,MVector(false),MVector(0.0),p, Dict("V"=>length(phase.species)+1)), y0, p
end
export ConstantTPDomain

struct ConstantVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    V::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantVDomain(;phase::Z,initialconds::Dict{X,E},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {E,X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert V != 0.0 || (T != 0.0 && P != 0.0)
    N = sum(ns)
    if V == 0.0
        V = N*R*T/P
    elseif T == 0.0
        T = P*V/(R*N)
    elseif P == 0.0
        P = N*R*T/V
    else
        throw(error("ConstantVDomain overspecified with T,P and V"))
    end
    y0 = vcat(ns,T,P)
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    else
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    end
    rxnarray = getreactionindices(phase)
    return ConstantVDomain(phase,[1,length(phase.species),length(phase.species)+1,length(phase.species)+2],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    V,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict("T"=>length(phase.species)+1,"P"=>length(phase.species)+2)), y0, p
end
export ConstantVDomain

struct ConstantPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    P::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantPDomain(;phase::Z,initialconds::Dict{X,E},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {E,X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert P != 0.0 || (T != 0.0 && V != 0.0)
    N = sum(ns)
    if P == 0.0
        P = N*R*T/V
    elseif T == 0.0
        T = P*V/(R*N)
    elseif V == 0.0
        V = N*R*T/P
    else
        throw(error("ConstantPDomain overspecified with T,P and V"))
    end
    y0 = vcat(ns,T,V)
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    else
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    end
    rxnarray = getreactionindices(phase)
    return ConstantPDomain(phase,[1,length(phase.species),length(phase.species)+1,length(phase.species)+2],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    P,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict("T"=>length(phase.species)+1,"V"=>length(phase.species)+2)), y0, p
end
export ConstantPDomain

struct ParametrizedTPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::Function
    P::Function
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedTPDomain(;phase::Z,initialconds::Dict{X,Any},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {X,X2,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = 0.0
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind) <: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert V != 0.0 || (T != 0.0 && P != 0.0)
    if isa(T,AbstractArray)
        Tfcn = getspline(ts,T)
    elseif isa(T,Function)
        Tfcn = T
    else
        throw(error("ParametrizedTPDomain must take \"T\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    if isa(P,AbstractArray)
        Pfcn = getspline(ts,P)
    elseif isa(P,Function)
        Pfcn = P
    else
        throw(error("ParametrizedTPDomain must take \"P\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end

    N = sum(ns)
    V = N*R*Tfcn(0.0)/Pfcn(0.0)
    y0 = zeros(length(phase.species)+1)
    y0[1:length(phase.species)] = ns
    y0[length(phase.species)+1] = V
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(V),length(phase.species)+1,length(phase.species)+1)
    else
        jacobian=zeros(typeof(V),length(phase.species)+1,length(phase.species)+1)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedTPDomain(phase,[1,length(phase.species),length(phase.species)+1],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    Tfcn,Pfcn,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict("V"=>length(phase.species)+1)), y0, p
end
export ParametrizedTPDomain

struct ParametrizedVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    V::Function
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedVDomain(;phase::Z,initialconds::Dict{X,Any},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {X,X2,E<:Real,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "V" in keys(initialconds)
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert isa(V,Function) || isa(V,AbstractArray)
    if isa(V,AbstractArray)
        Vfcn = getspline(ts,V)
    elseif isa(V,Function)
        Vfcn = V
    else
        throw(error("ParametrizedVDomain must take \"V\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    if T == 0.0
        T = P*Vfcn(0.0)/(R*N)
    elseif P == 0.0
        P = N*R*T/Vfcn(0.0)
    else
        ns *= (P*Vfcn(0.0)/(R*T))/sum(ns) #automatically scale down moles if pressure specified
    end
    y0 = vcat(ns,T,P)
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    else
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedVDomain(phase,[1,length(phase.species),length(phase.species)+1,length(phase.species)+2],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    Vfcn,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict("T"=>length(phase.species)+1,"P"=>length(phase.species)+2)), y0, p
end
export ParametrizedVDomain

struct ParametrizedPDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    P::Function
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedPDomain(;phase::Z,initialconds::Dict{X,Any},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {X,X2,E<:Real,Z<:IdealGas,Q<:AbstractInterface}

    #set conditions and initialconditions
    T = 0.0
    P = 0.0
    V = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "P" in keys(initialconds)
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            P = val
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    @assert isa(P,Function) || isa(P,AbstractArray)
    if isa(P,AbstractArray)
        Pfcn = getspline(ts,P)
    elseif isa(P,Function)
        Pfcn = P
    else
        throw(error("ParametrizedPDomain must take \"P\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    if T == 0.0
        T = Pfcn(0.0)*V/(R*N)
    elseif V == 0.0
        V = N*R*T/Pfcn(0.0)
    else
        ns *= (Pfcn(0.0)*V/(R*T))/sum(ns) #automatically scale down moles if volume specified
    end
    y0 = vcat(ns,T,V)
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    else
        jacobian=zeros(typeof(T),length(phase.species)+2,length(phase.species)+2)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedPDomain(phase,[1,length(phase.species),length(phase.species)+1,length(phase.species)+2],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    Pfcn,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict("T"=>length(phase.species)+1,"V"=>length(phase.species)+2)), y0, p
end
export ParametrizedPDomain

mutable struct ConstantTVDomain{N<:AbstractPhase,S<:Integer,W<:Real, W2<:Real, I<:Integer, Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    V::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    kfsnondiff::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantTVDomain(;phase::Z,initialconds::Dict{X,E},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse=false,sensitivity=false) where {E,X,X2, Z<:AbstractPhase,Q<:AbstractInterface,W<:Real}
    #set conditions and initialconditions
    T = 0.0
    V = 0.0
    P = 1.0e8
    y0 = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            throw(error("ConstantTVDomain cannot specify P"))
        elseif key == "V"
            V = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end
    @assert T != 0.0
    @assert V != 0.0
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase,T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    if phase.diffusionlimited
        diffs = [x(T=T,mu=mu,P=P) for x in getfield.(phase.species,:diffusion)]
    else
        diffs = Array{Float64,1}()
    end
    P = 1.0e8  #essentiallly assuming this is a liquid
    C = N/V
    kfs,krevs = getkfkrevs(phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    kfsnondiff = getkfs(phase,T,P,C,ns,V,0.0)
    p = vcat(Gs,kfsnondiff)
    if sparse
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    else
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTVDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
        T,V,kfs,krevs,kfsnondiff,efficiencyinds,Gs,rxnarray,mu,diffs,jacobian,sensitivity,false,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ConstantTVDomain

struct ParametrizedTConstantVDomain{N<:AbstractPhase,S<:Integer,W<:Real,W2<:Real,I<:Integer,Q<:AbstractArray} <: AbstractVariableKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::Function
    V::W
    efficiencyinds::Array{I,1}
    rxnarray::Array{Int64,2}
    jacobian::Array{W,2}
    sensitivity::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ParametrizedTConstantVDomain(;phase::IdealDiluteSolution,initialconds::Dict{X,X3},constantspecies::Array{X2,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false) where {X,X2,X3,Q<:AbstractInterface}
    #set conditions and initialconditions
    T = 0.0
    P = 1.0e8 #essentiallly assuming this is a liquid
    V = 0.0
    ts = Array{Float64,1}()
    ns = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    @assert "V" in keys(initialconds)
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "P"
            throw(error("ParametrizedTConstantVDomain cannot specify P"))
        elseif key == "V"
            V = val
        elseif key == "ts"
            ts = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            ns[ind] = val
        end
    end
    if isa(T,AbstractArray)
        Tfcn = getspline(ts,T)
    elseif isa(T,Function)
        Tfcn = T
    else
        throw(error("ParametrizedTConstantVDomain must take \"T\" as a function or if an array of times for \"ts\" is supplied as an array of volumes"))
    end
    N = sum(ns)
    y0 = zeros(length(phase.species))
    y0[1:length(phase.species)] = ns
    p = vcat(zeros(length(phase.species)),ones(length(phase.reactions)))
    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    if sparse
        jacobian=zeros(typeof(V),length(phase.species)+1,length(phase.species)+1)
    else
        jacobian=zeros(typeof(V),length(phase.species)+1,length(phase.species)+1)
    end
    rxnarray = getreactionindices(phase)
    return ParametrizedTConstantVDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
    Tfcn,V,efficiencyinds,rxnarray,jacobian,sensitivity,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ParametrizedTConstantVDomain

mutable struct ConstantTAPhiDomain{N<:AbstractPhase,S<:Integer,W<:Real, W2<:Real, I<:Integer, Q<:AbstractArray} <: AbstractConstantKDomain
    phase::N
    indexes::Q #assumed to be in ascending order
    parameterindexes::Q
    constantspeciesinds::Array{S,1}
    T::W
    A::W
    phi::W
    kfs::Array{W,1}
    krevs::Array{W,1}
    efficiencyinds::Array{I,1}
    Gs::Array{W,1}
    rxnarray::Array{Int64,2}
    mu::W
    diffusivity::Array{W,1}
    jacobian::Array{W,2}
    sensitivity::Bool
    alternativepformat::Bool
    jacuptodate::MArray{Tuple{1},Bool,1,1}
    t::MArray{Tuple{1},W2,1,1}
    p::Array{W,1}
    thermovariabledict::Dict{String,Int64}
end
function ConstantTAPhiDomain(;phase::E2,initialconds::Dict{X,X2},constantspecies::Array{X3,1}=Array{String,1}(),
    sparse::Bool=false,sensitivity::Bool=false,stationary::Bool=false) where {E<:Real,E2<:AbstractPhase,W<:Real,X,X2,X3}
    #set conditions and initialconditions
    T = 0.0
    A = 0.0
    phi = 0.0 #default 0.0
    y0 = zeros(length(phase.species))
    spnames = [x.name for x in phase.species]
    for (key,val) in initialconds
        if key == "T"
            T = val
        elseif key == "A"
            A = val
        elseif key == "Phi"
            phi = val
        else
            ind = findfirst(isequal(key),spnames)
            @assert typeof(ind)<: Integer  "$key not found in species list: $spnames"
            y0[ind] = val
        end
    end

    @assert A != 0.0
    @assert T != 0.0
    ns = y0
    N = sum(ns)

    if length(constantspecies) > 0
        spcnames = getfield.(phase.species,:name)
        constspcinds = [findfirst(isequal(k),spcnames) for k in constantspecies]
    else
        constspcinds = Array{Int64,1}()
    end
    efficiencyinds = [rxn.index for rxn in phase.reactions if typeof(rxn.kinetics)<:AbstractFalloffRate && length(rxn.kinetics.efficiencies) > 0]
    Gs = calcgibbs(phase,T)
    if :solvent in fieldnames(typeof(phase)) && typeof(phase.solvent) != EmptySolvent
        mu = phase.solvent.mu(T)
    else
        mu = 0.0
    end
    C = 0.0 #this currently shouldn't matter here, on a surface you shouldn't have pdep
    kfs,krevs = getkfkrevs(phase,T,0.0,C,N,ns,Gs,[],A,phi)
    p = vcat(Gs,kfs)
    if sparse
        jacobian=spzeros(typeof(T),length(phase.species),length(phase.species))
    else
        jacobian=zeros(typeof(T),length(phase.species),length(phase.species))
    end
    rxnarray = getreactionindices(phase)
    return ConstantTAPhiDomain(phase,[1,length(phase.species)],[1,length(phase.species)+length(phase.reactions)],constspcinds,
        T,A,phi,kfs,krevs,efficiencyinds,Gs,rxnarray,mu,Array{Float64,1}(),jacobian,sensitivity,false,MVector(false),MVector(0.0),p,Dict{String,Int64}()), y0, p
end
export ConstantTAPhiDomain

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q,p::W3=DiffEqBase.NullParameters()) where {W3<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:Array{Float64,1},Q} #no parameter input
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P*V/(R*d.T)
    cs = ns./V
    C = N/V
    for ind in d.efficiencyinds #efficiency related rates may have changed
        d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0)
    end
    return ns,cs,d.T,d.P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:Array{Float64,1},Q<:Float64} #uses parameter input
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P*V/(R*d.T)
    cs = ns./V
    C = N/V
    if !d.alternativepformat
        @views kfps = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        @views nothermochg= d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else 
        kfps = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        nothermochg= d.Gs == d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    @views nokfchg = count(d.kfs .!= kfps) <= length(d.efficiencyinds) && all(kfps[d.efficiencyinds] .== 1.0)
    if nothermochg && nokfchg
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0;f=kfps[ind])
        end
        return ns,cs,d.T,d.P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
    elseif nothermochg
        d.kfs = kfps
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0;f=kfps[ind]) 
        end
        return ns,cs,d.T,d.P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
    else #need to handle thermo changes
        d.kfs .= kfps
        if !d.alternativepformat
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        else
            d.Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        end
        krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0;kfs=d.kfs)[2]
        for ind in d.efficiencyinds #efficiency related rates may have changed
            d.kfs[ind],d.krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,d.Gs,d.diffusivity,V,0.0;f=kfps[ind])
        end
        return ns,cs,d.T,d.P,V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
    end
end

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::Array{W3,1},t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,W3<:ForwardDiff.Dual,Q} #Autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P*V/(R*d.T)
    cs = ns./V
    C = N/V
    if !d.alternativepformat
        kfs = convert(typeof(y),p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else 
        kfs = convert(typeof(y),d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    krevs = convert(typeof(y),getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0;kfs=kfs)[2])
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0;f=kfs[ind])
    end
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J,Q} #Autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P*V/(R*d.T)
    cs = ns./V
    C = N/V
    if !d.alternativepformat
        kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    else 
        kfs = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
    end
    krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0;kfs=kfs)[2]
    for ind in d.efficiencyinds #efficiency related rates may have changed
        kfs[ind],krevs[ind] = getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0;f=kfs[ind])
    end
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Union{ReverseDiff.TrackedArray,Tracker.TrackedArray},W<:IdealGas,Y<:Integer,J,Q} #Tracker/reversediff
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = d.P*V/(R*d.T)
    cs = ns./V
    C = N/V
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0)[1]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    else 
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = [ind in d.efficiencyinds ? getkfkrev(d.phase.reactions[ind],d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0)[1]*d.p[length(d.phase.species)+ind]*p[d.parameterindexes[1]-1+length(d.phase.species)+ind] : p[d.parameterindexes[1]-1+length(d.phase.species)+ind] for ind in 1:length(d.phase.reactions)]
    end
    krevs = getkfkrevs(d.phase,d.T,d.P,C,N,ns,Gs,d.diffusivity,V,0.0;kfs=kfs)[2]
    return ns,cs,d.T,d.P,V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*d.V/(R*T)
    cs = ns./d.V
    C = N/d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0)
    return ns,cs,T,P,d.V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*d.V/(R*T)
    cs = ns./d.V
    C = N/d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @views @fastmath hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0)
    return @views @fastmath ns,cs,T,P,d.V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*d.V/(R*T)
    cs = ns./d.V
    C = N/d.V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    Cvave = 0.0
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,d.V,0.0)
    return @views @fastmath ns,cs,T,P,d.V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=d.P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0)
    return ns,cs,T,d.P,V,C,N,0.0,kfs,krevs,Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ConstantPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=d.P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0)
    if p != DiffEqBase.NullParameters()
        return @views @fastmath ns,cs,T,d.P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
    else
        return ns,cs,T,d.P,V,C,N,0.0,kfs,krevs,Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
    end
end

@inline function calcthermo(d::ConstantPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = d.P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=d.P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,d.P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,d.P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    V = d.V(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    P = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Us = (hdivRT.-1.0)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Us,Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end
@inline function calcthermo(d::ParametrizedPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    T = y[d.indexes[3]]
    V = y[d.indexes[4]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Hs = zeros(length(d.phase.species))
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Hs = hdivRT.*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Hs,Array{Float64,1}(),Gs,diffs,Cvave,cpdivR,0.0
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./V
    C = N/V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    if d.phase.diffusionlimited
        diffs = [x(T=T,mu=mu,P=P) for x in getfield.(d.phase.species,:diffusion)]::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return ns,cs,T,P,V,C,N,mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./V
    C = N/V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    if d.phase.diffusionlimited
        diffs = [x(T=T,mu=mu,P=P) for x in getfield.(d.phase.species,:diffusion)]::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,mu,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTConstantVDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    V = d.V
    T = d.T(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./V
    C = N/V
    P = 1.0e8 #liquid phase
    Gs = zeros(length(d.phase.species))
    mu = d.phase.solvent.mu(T)
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    if d.phase.diffusionlimited
        diffs = [x(T=T,mu=mu,P=P) for x in getfield.(d.phase.species,:diffusion)]::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,mu,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:DiffEqBase.NullParameters,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return ns,cs,T,P,V,C,N,0.0,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2<:Array{Float64,1},W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT .+= p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ParametrizedTPDomain{W,Y},y::J,t::Q,p::W2=DiffEqBase.NullParameters()) where {W2,W<:IdealGas,Y<:Integer,J<:AbstractArray,Q}
    T = d.T(t)
    @assert T < 10000.0
    P = d.P(t)
    ns = y[d.indexes[1]:d.indexes[2]]
    V = y[d.indexes[3]]
    N = P*V/(R*T)
    cs = ns./V
    C = N/V
    Gs = zeros(length(d.phase.species))
    Us = zeros(length(d.phase.species))
    cpdivR,hdivRT1,sdivR = calcHSCpdless(d.phase.vecthermo,T)
    @fastmath @views hdivRT = hdivRT1 .+ p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]./(R*T)
    @fastmath Gs = (hdivRT.-sdivR)*(R*T)
    @fastmath Cvave = dot(cpdivR,ns)
    @fastmath Cvave *= R/N
    @fastmath Cvave -= R
    if d.phase.diffusionlimited
        diffs = getfield.(d.phase.species,:diffusion)(T=T,mu=0.0,P=P)::Array{typeof(T),1}
    else
        diffs = Array{Float64,1}()
    end
    kfs,krevs = getkfkrevs(d.phase,T,P,C,N,ns,Gs,diffs,V,0.0)
    return @views @fastmath ns,cs,T,P,V,C,N,0.0,kfs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],krevs.*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(kfs)],Array{Float64,1}(),Array{Float64,1}(),Gs,diffs,0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTVDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2<:DiffEqBase.NullParameters,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = 1.0e8
    return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTVDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2<:Array{Float64,1},W<:IdealDiluteSolution,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = 1.0e8
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        elseif nothermochg
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        else
            d.kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        end
    else 
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        @views nokfchg = d.kfsnondiff == d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
        if nothermochg && nokfchg
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        elseif nothermochg
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        else
            d.kfsnondiff .= d.p[length(d.phase.species)+1:length(d.phase.species)+length(d.phase.reactions)].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs .= d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.kfs,d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,0.0;kfs=d.kfsnondiff)
            return ns,cs,d.T,P,d.V,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
        end
    end
        
end

@inline function calcthermo(d::ConstantTVDomain{W,Y},y::Array{W2,1},t::Q,p::Q2=DiffEqBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y),p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else 
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = convert(typeof(y),d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.V,0.0;kfs=kfsnondiff)
    return ns,cs,d.T,P,d.V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTVDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2,W<:IdealDiluteSolution,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.V
    C = N/d.V
    P = 1.0e8
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else 
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfsnondiff = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]    
    end
    kfs,krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.V,0.0;kfs=kfsnondiff)
    return ns,cs,d.T,P,d.V,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),Array{Float64,1}(),0.0,Array{Float64,1}(),0.0
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2<:DiffEqBase.NullParameters,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.A
    C = N/d.A
    P = 0.0
    return ns,cs,d.T,P,d.A,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2<:Array{Float64,1},W<:IdealSurface,Y<:Integer,J<:Array{Float64,1},Q}
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.A
    C = N/d.A
    P = 0.0
    if !d.alternativepformat
        @views nothermochg = d.Gs == p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        if nothermochg
            return ns,cs,d.T,P,d.A,C,N,d.mu,p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)],d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        else
            d.kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi;kfs=d.kfs)[2]
            return ns,cs,d.T,P,d.A,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    else 
        @views nothermochg = d.Gs == d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        if nothermochg
            return ns,cs,d.T,P,d.A,C,N,d.mu,d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)],d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,0.0,Array{Float64,1}(),d.phi
        else
            d.kfs = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
            d.Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
            d.krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,d.Gs,d.diffusivity,d.V,d.phi;kfs=d.kfs)[2]
            return ns,cs,d.T,P,d.A,C,N,d.mu,d.kfs,d.krevs,Array{Float64,1}(),Array{Float64,1}(),d.Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
        end
    end
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y},y::Array{W2,1},t::Q,p::Q2=DiffEqBase.NullParameters()) where {W2<:ForwardDiff.Dual,Q2,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q} #autodiff y
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.A
    C = N/d.A
    P = 0.0
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = convert(typeof(y),p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    else 
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = convert(typeof(y),d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)])
    end
    krevs = convert(typeof(y),getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.A,d.phi;kfs=kfs)[2])
    return ns,cs,d.T,P,d.A,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end

@inline function calcthermo(d::ConstantTAPhiDomain{W,Y},y::J,t::Q,p::Q2=DiffEqBase.NullParameters()) where {Q2,W<:IdealSurface,Y<:Integer,J<:AbstractArray,Q} #autodiff p
    ns = y[d.indexes[1]:d.indexes[2]]
    N = sum(ns)
    cs = ns./d.A
    C = N/d.A
    P = 0.0
    if !d.alternativepformat
        Gs = p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    else 
        Gs = d.p[1:length(d.phase.species)].+p[d.parameterindexes[1]-1+1:d.parameterindexes[1]-1+length(d.phase.species)]
        kfs = d.p[length(d.phase.species)+1:end].*p[d.parameterindexes[1]-1+length(d.phase.species)+1:d.parameterindexes[1]-1+length(d.phase.species)+length(d.phase.reactions)]
    end
    krevs = getkfkrevs(d.phase,d.T,P,C,N,ns,Gs,d.diffusivity,d.A,d.phi;kfs=kfs)[2]
    return ns,cs,d.T,P,d.A,C,N,d.mu,kfs,krevs,Array{Float64,1}(),Array{Float64,1}(),Gs,Array{Float64,1}(),0.0,Array{Float64,1}(),d.phi
end
export calcthermo

@inline function calcdomainderivatives!(d::Q,dydt::Z7,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Array{Z,1},Hs::Array{Z11,1},V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z8) where {Q<:AbstractDomain,Z12,Z11,Z10,Z9,Z8<:Real,Z7,W<:IdealGas,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*inter.F(t)
        elseif isa(inter,Outlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.F(t).*ns./N
        elseif isa(inter,Condensation) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.kLA./inter.kH.*inter.x_vap*inter.P_vap*V/R/inter.T_vap
        elseif isa(inter,Evaporation) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.kLA./inter.kH.*ns
        end
    end
end

@inline function calcdomainderivatives!(d::Q,dydt::Z7,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Array{Z,1},Hs::Array{Z11,1},V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z8) where {Q<:ConstantTPDomain,Z12,Z11,Z10,Z9,Z8<:Real,Z7,W<:IdealGas,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @views @fastmath @inbounds dydt[d.indexes[3]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/P
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*inter.F(t)
            dydt[d.indexes[3]] += inter.F(t)*R*T/P
        elseif isa(inter,Outlet) && d == inter.domain
            dydt[d.indexes[1]:d.indexes[2]] .-= inter.F(t).*ns./N
            dydt[d.indexes[3]] -= inter.F(t)*R*T/P
        end
    end
end

@inline function calcdomainderivatives!(d::ConstantVDomain{W,Y},dydt::K,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Z,Hs::Z11,V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z7) where {Z12,Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @views @fastmath @inbounds dydt[d.indexes[3]] = -dot(Us,dydt[d.indexes[1]:d.indexes[2]])/(N*Cvave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/V + P/T*dydt[d.indexes[3]]
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*flow
            dTdt = flow*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow*R*T/V + P/T*dTdt
        elseif isa(inter,Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow.*ns./N
            dTdt = (P*V/N*flow)/(N*Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow*R*T/V + P/T*dTdt
        end
    end
end

@inline function calcdomainderivatives!(d::ConstantPDomain{W,Y},dydt::K,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Z,Hs::Z11,V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z7) where {Z12,Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5}
    @fastmath Cpave = Cvave+R
    @views @fastmath @inbounds dydt[d.indexes[3]] = -dot(Hs,dydt[d.indexes[1]:d.indexes[2]])/(N*Cpave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/P + dydt[d.indexes[3]]*V/T
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*flow
            dTdt = flow*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow*R*T/P + dTdt*V/T 
        elseif isa(inter,Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow.*ns./N
            dydt[d.indexes[4]] -= flow*R*T/P
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedTPDomain{W,Y},dydt::K,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Z,Hs::Z11,V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    @views @fastmath @inbounds dydt[d.indexes[3]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/P + Calculus.derivative(d.T,t)*V/T - Calculus.derivative(d.P,t)*V/P
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*flow
            dydt[d.indexes[3]] += flow*R*T/P
        elseif isa(inter,Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow*ns./N
            dydt[d.indexes[3]] -= flow*R*T/P
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedVDomain{W,Y},dydt::K,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Z,Hs::Z11,V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    dVdt = Calculus.derivative(d.V,t)
    @views @fastmath @inbounds dydt[d.indexes[3]] = (-dot(Us,dydt[d.indexes[1]:d.indexes[2]])-P*dVdt)/(N*Cvave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/V + dydt[d.indexes[3]]*P/T - P/V*dVdt
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*flow
            dTdt = flow*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow*R*T/V + dTdt*P/T
        elseif isa(inter,Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow*ns./N
            dTdt = (P*V/N*flow)/(N*Cvave)
            dydt[d.indexes[3]] -= dTdt
            dydt[d.indexes[4]] -= flow*R*T/V + dTdt*P/T
        end
    end
end

@inline function calcdomainderivatives!(d::ParametrizedPDomain{W,Y},dydt::K,interfaces::Z12;t::Z10,T::Z4,P::Z9,Us::Z,Hs::Z11,V::Z2,C::Z3,ns::Z5,N::Z6,Cvave::Z7) where {Z11,Z10,Z9,W<:IdealGas,Z7,K,Y<:Integer,Z6,Z,Z2,Z3,Z4,Z5,Z12}
    @fastmath Cpave = Cvave+R
    dPdt = Calculus.derivative(d.P,t)
    @views @fastmath @inbounds dydt[d.indexes[3]] = (-dot(Hs,dydt[d.indexes[1]:d.indexes[2]])+V*dPdt)/(N*Cpave) #divide by V to cancel ωV to ω
    @views @fastmath @inbounds dydt[d.indexes[4]] = sum(dydt[d.indexes[1]:d.indexes[2]])*R*T/P + dydt[d.indexes[3]]*V/T - dPdt*V/P 
    for ind in d.constantspeciesinds #make dydt zero for constant species
        @inbounds dydt[ind] = 0.0
    end
    for inter in interfaces
        if isa(inter,Inlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .+= inter.y.*flow
            dTdt = flow*(inter.H - dot(Hs,ns)/N)/(N*Cpave)
            dydt[d.indexes[3]] += dTdt
            dydt[d.indexes[4]] += flow*R*T/P + dTdt*V/T
        elseif isa(inter,Outlet) && d == inter.domain
            flow = inter.F(t)
            dydt[d.indexes[1]:d.indexes[2]] .-= flow.*ns./N
            dydt[d.indexes[4]] -= flow*R*T/P
        end
    end
end
export calcdomainderivatives!

@inline function jacobianyefficiencyderiv!(jac::S,domain::Union{ConstantTPDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedPDomain},kinetics::AbstractFalloffRate,efficiencies::Dict{Int64,Float64},rxnarray::Array{Int64,2},rxnind::Int64,Vind::Int64,cs::Array{Float64,1},T::Float64,V::Float64,C::Float64,Ceff::Float64,Kc::Float64) where {S<:AbstractArray}
    dkdCeff = _calcdkdCeff(kinetics,T,Ceff)
    for (spcind,efficiencyval) in efficiencies
        @fastmath dkdni = dkdCeff*(efficiencyval/V)
        _jacobianykderiv!(jac,spcind,dkdni,dkdni/Kc,rxnarray,rxnind,cs,V)
    end
    @fastmath dkdV = dkdCeff*(C-Ceff)/V
    _jacobianykderiv!(jac,Vind,dkdV,dkdV/Kc,rxnarray,rxnind,cs,V)
end

@inline function jacobianyefficiencyderiv!(jac::S,domain::Union{ConstantVDomain,ParametrizedVDomain,ConstantTVDomain,ParametrizedTConstantVDomain},kinetics::AbstractFalloffRate,efficiencies::Dict{Int64,Float64},rxnarray::AbstractArray,rxnind::Int64,Vind::Int64,cs::Array{Float64,1},T::Float64,V::Float64,C::Float64,Ceff::Float64,Kc::Float64) where {S<:AbstractArray}
    dkdCeff = _calcdkdCeff(kinetics,T,Ceff)
    for (spcind,efficiencyval) in efficiencies
        @fastmath dkdni = dkdCeff*(efficiencyval/V)
        _jacobianykderiv!(jac,spcind,dkdni,dkdni/Kc,rxnarray,rxnind,cs,V)
    end
end

@inline function jacobianynsderiv!(jac::S,domain::Union{ConstantTPDomain,ParametrizedTPDomain},rxnarray::Array{Int64,2},efficiencyinds::Array{I,1},cs::Array{Float64,1},kfs::Array{Float64,1},krevs::Array{Float64,1},T::Float64,V::Float64,C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac,rxnarray,rxnind,cs,kfs[rxnind],krevs[rxnind])
        @inbounds _jacobianynswrtV!(jac,domain.indexes[3],rxnarray,rxnind,cs,kfs[rxnind],krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind]/krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i]*val for (i,val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac,domain,domain.phase.reactions[rxnind].kinetics,domain.phase.reactions[rxnind].kinetics.efficiencies,rxnarray,rxnind,domain.indexes[3],cs,T,V,C,Ceff,Kc)
        end
    end
end

@inline function jacobianynsderiv!(jac::S,domain::Union{ConstantVDomain,ParametrizedVDomain,ConstantTVDomain,ConstantTAPhiDomain,ParametrizedTConstantVDomain},rxnarray::Array{Int64,2},efficiencyinds::Array{I,1},cs::Array{Float64,1},kfs::Array{Float64,1},krevs::Array{Float64,1},T::Float64,V::Float64,C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac,rxnarray,rxnind,cs,kfs[rxnind],krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind]/krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i]*val for (i,val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac,domain,domain.phase.reactions[rxnind].kinetics,domain.phase.reactions[rxnind].kinetics.efficiencies,rxnarray,rxnind,0,cs,T,V,C,Ceff,Kc)
        end
    end
end

@inline function jacobianynsderiv!(jac::S,domain::Union{ConstantPDomain,ParametrizedPDomain},rxnarray::Array{Int64,2},efficiencyinds::Array{I,1},cs::Array{Float64,1},kfs::Array{Float64,1},krevs::Array{Float64,1},T::Float64,V::Float64,C::Float64) where {S<:AbstractArray,I<:Integer}
    @simd for rxnind = 1:size(rxnarray)[2]
        @inbounds _jacobianynswrtns!(jac,rxnarray,rxnind,cs,kfs[rxnind],krevs[rxnind])
        @inbounds _jacobianynswrtV!(jac,domain.indexes[4],rxnarray,rxnind,cs,kfs[rxnind],krevs[rxnind])
        if rxnind in efficiencyinds
            @inbounds @fastmath Kc = kfs[rxnind]/krevs[rxnind]
            @inbounds @fastmath Ceff = C + sum(cs[i]*val for (i,val) in domain.phase.reactions[rxnind].kinetics.efficiencies)
            @inbounds jacobianyefficiencyderiv!(jac,domain,domain.phase.reactions[rxnind].kinetics,domain.phase.reactions[rxnind].kinetics.efficiencies,rxnarray,rxnind,domain.indexes[4],cs,T,V,C,Ceff,Kc)
        end
    end
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantTPDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @views @inbounds @fastmath jac[domain.indexes[3],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/P
    end
    @views @inbounds @fastmath jac[domain.indexes[3],domain.indexes[3]] = sum(jac[domain.indexes[1]:domain.indexes[2],domain.indexes[3]])*R*T/P
          
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2],domain.indexes[3]] .-= flow.*ns.*(-C)
        end
    end
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantVDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)

    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V

    for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCvavedni = cpdivR[i]*R/N
        @views @inbounds @fastmath jac[domain.indexes[3],i] = -dot(Us,jac[domain.indexes[1]:domain.indexes[2],i])/(N*Cvave)+dot(Us,dydt[domain.indexes[1]:domain.indexes[2]])/(N*Cvave*Cvave)*dCvavedni
        @views @inbounds @fastmath jac[domain.indexes[4],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/V + P/T*jac[domain.indexes[3],i]
    end
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i]*R/N
                @fastmath dTdt = flow*(inter.H - dot(Us,ns)/N)/(N*Cvave)
                @inbounds @fastmath ddnidTdt = flow*(-Us[i]/N)/(N*Cvave)-dTdt*(dCvavedni/Cvave)
                @inbounds jac[domain.indexes[3],i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] += P/T*ddnidTdt
            end
        elseif isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
                @inbounds @fastmath dCvavedni = cpdivR[i]*R/N
                @fastmath dTdt = (P*V/N*flow)/(N*Cvave)
                @fastmath ddnidTdt = -dTdt*(dCvavedni/Cvave)
                @inbounds jac[domain.indexes[3],i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] -= P/T*ddnidTdt
            end
        end
    end
    
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[3],T,colorvec)
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[4],P,colorvec)

    return jac
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ConstantPDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V
    
    @fastmath Cpave = Cvave+R
    @views @inbounds @fastmath dTdt = -dot(Hs,dydt[domain.indexes[1]:domain.indexes[2]])/(N*Cpave)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCpavedni = cpdivR[i]*R/N
        @views @inbounds @fastmath jac[domain.indexes[3],i] = -dot(Hs,jac[domain.indexes[1]:domain.indexes[2],i])/(N*Cpave)-dTdt*(dCpavedni/Cpave)
        @views @inbounds @fastmath jac[domain.indexes[4],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/P + V/T*jac[domain.indexes[3],i]
    end
    @views @inbounds @fastmath jac[domain.indexes[3],domain.indexes[4]] = -dot(Hs,jac[domain.indexes[1]:domain.indexes[2],domain.indexes[4]])/(N*Cpave)
    @views @inbounds @fastmath jac[domain.indexes[4],domain.indexes[4]] = sum(jac[domain.indexes[1]:domain.indexes[2],domain.indexes[4]])*R*T/P + dTdt/T + V/T*jac[domain.indexes[3],domain.indexes[4]]
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            @fastmath H = dot(Hs,ns)/N
            @fastmath dTdt = flow*(inter.H - H)/(N*Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i]*R/N
                @inbounds @fastmath ddnidTdt = flow*(-Hs[i]/N)/(N*Cpave)-dTdt*(dCpavedni/Cpave)
                @inbounds jac[domain.indexes[3],i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] += V/T*ddnidTdt
            end
            @fastmath ddVdTdt = -flow*H/V/(N*Cpave)
            @inbounds jac[domain.indexes[3],domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4],domain.indexes[4]] += dTdt/T + V/T*ddVdTdt
        elseif isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
            end
        end
    end
    
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[3],T,colorvec)
    
    return jac
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedTPDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @views @inbounds @fastmath jac[domain.indexes[3],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/P
    end
    @views @inbounds @fastmath jac[domain.indexes[3],domain.indexes[3]] = sum(jac[domain.indexes[1]:domain.indexes[2],domain.indexes[3]])*R*T/P + Calculus.derivative(domain.T,t)/T - Calculus.derivative(domain.P,t)/P
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
            end
        end
    end
    
    return jac
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedVDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V
    
    dVdt = Calculus.derivative(domain.V,t)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @fastmath @inbounds dCvavedni = cpdivR[i]*R/N
        @views @inbounds @fastmath jac[domain.indexes[3],i] = -dot(Us,jac[domain.indexes[1]:domain.indexes[2],i])/(N*Cvave)-(-dot(Us,dydt[domain.indexes[1]:domain.indexes[2]])-P*dVdt)/(N*Cvave*Cvave)*dCvavedni
        @views @inbounds @fastmath jac[domain.indexes[4],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/V + P/T*jac[domain.indexes[3],i]
    end
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            @fastmath dTdt = flow*(inter.H - dot(Us,ns)/N)/(N*Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCvavedni = cpdivR[i]*R/N
                @inbounds @fastmath ddnidTdt = flow*(-Us[i]/N)/(N*Cvave)-dTdt*(dCvavedni/Cvave)
                @inbounds jac[domain.indexes[3],i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] += P/T*ddnidTdt
            end
        elseif isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @fastmath dTdt = (P*V/N*flow)/(N*Cvave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
                @inbounds @fastmath dCvavedni = cpdivR[i]*R/N
                @fastmath ddnidTdt = -dTdt*(dCvavedni/Cvave)
                @inbounds jac[domain.indexes[3],i] -= ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] -= P/T*ddnidTdt
            end
        end
    end
    
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[3],T,colorvec)
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[4],P,colorvec)

    return jac
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:ParametrizedPDomain}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V
    
    @fastmath Cpave = Cvave+R
    dPdt = Calculus.derivative(domain.P,t)
    @views @inbounds @fastmath dTdt = (-dot(Hs,dydt[domain.indexes[1]:domain.indexes[2]])+V*dPdt)/(N*Cpave)
    @simd for i in domain.indexes[1]:domain.indexes[2]
        @inbounds @fastmath dCpavedni = cpdivR[i]*R/N
        @views @inbounds @fastmath jac[domain.indexes[3],i] = -dot(Hs,jac[domain.indexes[1]:domain.indexes[2],i])/(N*Cpave)-dTdt*(dCpavedni/Cpave)
        @views @inbounds @fastmath jac[domain.indexes[4],i] = sum(jac[domain.indexes[1]:domain.indexes[2],i])*R*T/P + V/T*jac[domain.indexes[3],i]
    end
    @views @inbounds @fastmath jac[domain.indexes[3],domain.indexes[4]] = (-dot(Hs,jac[domain.indexes[1]:domain.indexes[2],domain.indexes[4]])+dPdt)/(N*Cpave)
    @views @inbounds @fastmath jac[domain.indexes[4],domain.indexes[4]] = sum(jac[domain.indexes[1]:domain.indexes[2],domain.indexes[4]])*R*T/P + dTdt/T + V/T*jac[domain.indexes[3],domain.indexes[4]]-dPdt/P 
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            @fastmath H = dot(Hs,ns)/N
            @fastmath dTdt = flow*(inter.H - H)/(N*Cpave)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath dCpavedni = cpdivR[i]*R/N
                @inbounds @fastmath ddnidTdt = flow*(-Hs[i]/N)/(N*Cpave)-dTdt*(dCpavedni/Cpave)
                @inbounds jac[domain.indexes[3],i] += ddnidTdt
                @inbounds @fastmath jac[domain.indexes[4],i] += V/T*ddnidTdt
            end
            @fastmath ddVdTdt = -flow*H/V/(N*Cpave)
            @inbounds jac[domain.indexes[3],domain.indexes[4]] += ddVdTdt
            @inbounds @fastmath jac[domain.indexes[4],domain.indexes[4]] += dTdt/T + ddVdTdt*V/T
        elseif isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
            end
            @views @inbounds @fastmath jac[domain.indexes[1]:domain.indexes[2],domain.indexes[4]] .-= -flow.*ns./N/V
        end
    end
    
    @inbounds jacobianytherm!(jac,y,p,t,domain,interfaces,domain.indexes[3],T,colorvec)
    
    return jac
end

@inline function jacobiany!(jac::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain,ConstantTAPhiDomain,ParametrizedTConstantVDomain}}
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,phi = calcthermo(domain,y,t,p)
    jacobianynsderiv!(jac,domain,domain.rxnarray,domain.efficiencyinds,cs,kfs,krevs,T,V,C)
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jac[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Outlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= flow/N
            end
        elseif isa(inter,Evaporation) && domain == inter.domain
            @simd for i in domain.indexes[1]:domain.indexes[2]
                @inbounds @fastmath jac[i,i] -= inter.kLA[i]/inter.kH[i]
            end
        end
    end

    return jac
end

export jacobiany!

@inline function jacobianpnsderiv!(jacp::Q,y::U,p::W,t::Z,domain::D,rxnarray::Array{Int64,2},cs::Array{Float64,1},T::Float64,V::Float64,kfs::Array{Float64,1},krevs::Array{Float64,1},Nspcs::Int64,Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTPDomain,ConstantTAPhiDomain}}
    @fastmath  RTinv = 1.0/(R*T)
    @simd for rxnind = 1:Nrxns
        if @inbounds rxnarray[2,rxnind] == 0
            @inbounds fderiv = cs[rxnarray[1,rxnind]]
        elseif @inbounds  rxnarray[3,rxnind] == 0
            @fastmath @inbounds fderiv = cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
        else
            @fastmath @inbounds fderiv = cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]
        end

        if @inbounds rxnarray[5,rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind]/kfs[rxnind]*cs[rxnarray[4,rxnind]]
        elseif @inbounds rxnarray[6,rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind]/kfs[rxnind]*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
        else
            @fastmath @inbounds rderiv = krevs[rxnind]/kfs[rxnind]*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]
        end

        flux = fderiv-rderiv
        _spreadreactantpartials!(jacp,flux,rxnarray,rxnind,Nspcs+rxnind)
        _spreadproductpartials!(jacp,-flux,rxnarray,rxnind,Nspcs+rxnind)
        
        @fastmath @inbounds gderiv = rderiv*kfs[rxnind]*RTinv

        @inbounds jacp[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[1,rxnind])
        if @inbounds rxnarray[2,rxnind] !== 0
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[2,rxnind])
            if @inbounds rxnarray[3,rxnind] !== 0
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[1,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[2,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[3,rxnind])
            end
        end

        @inbounds jacp[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[4,rxnind])
        if @inbounds rxnarray[5,rxnind] !== 0
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[5,rxnind])
            if @inbounds rxnarray[6,rxnind] !== 0
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[4,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[5,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[4,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[6,rxnind])
            end
        end
        
    end
    jacp .*= V
end

@inline function jacobianpnsderiv!(jacp::Q,y::U,p::W,t::Z,domain::D,rxnarray::Array{Int64,2},cs::Array{Float64,1},T::Float64,V::Float64,kfs::Array{Float64,1},krevs::Array{Float64,1},Nspcs::Int64,Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantVDomain,ConstantPDomain,ParametrizedTPDomain,ParametrizedVDomain,ParametrizedPDomain,ParametrizedTConstantVDomain}}
    @fastmath  RTinv = 1.0/(R*T)
    @simd for rxnind = 1:Nrxns
        if @inbounds rxnarray[2,rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind]*cs[rxnarray[1,rxnind]]
        elseif @inbounds rxnarray[3,rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind]*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
        else
            @fastmath @inbounds fderiv = kfs[rxnind]*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]
        end

        if @inbounds rxnarray[5,rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind]*cs[rxnarray[4,rxnind]]
        elseif @inbounds rxnarray[6,rxnind] == 0
            @fastmath @inbounds rderiv = krevs[rxnind]*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
        else
            @fastmath @inbounds rderiv = krevs[rxnind]*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]
        end

        @fastmath flux = fderiv-rderiv
        _spreadreactantpartials!(jacp,flux,rxnarray,rxnind,Nspcs+rxnind)
        _spreadproductpartials!(jacp,-flux,rxnarray,rxnind,Nspcs+rxnind)
        
        @fastmath gderiv = rderiv*RTinv

        @inbounds jacp[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[1,rxnind])
        if @inbounds rxnarray[2,rxnind] !== 0
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[2,rxnind])
            if @inbounds rxnarray[3,rxnind] !== 0
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[1,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[2,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[3,rxnind])
            end
        end

        @inbounds jacp[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[4,rxnind])
        if @inbounds rxnarray[5,rxnind] !== 0
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[5,rxnind])
            if @inbounds rxnarray[6,rxnind] !== 0
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[4,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[5,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[4,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[6,rxnind])
            end
        end
    end
    jacp .*= V
end

@inline function jacobianpnsderiv!(jacp::Q,y::U,p::W,t::Z,domain::D,rxnarray::Array{Int64,2},cs::Array{Float64,1},T::Float64,V::Float64,kfs::Array{Float64,1},krevs::Array{Float64,1},Nspcs::Int64,Nrxns::Int64) where {Q3<:AbstractArray,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain}}
    @fastmath RTinv = 1.0/(R*T)
    @simd for rxnind = 1:Nrxns
        if @inbounds rxnarray[2,rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind]*kfs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[1,rxnind]]
        elseif @inbounds rxnarray[3,rxnind] == 0
            @fastmath @inbounds fderiv = kfs[rxnind]*kfs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]
        else
            @fastmath @inbounds fderiv = kfs[rxnind]*kfs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[1,rxnind]]*cs[rxnarray[2,rxnind]]*cs[rxnarray[3,rxnind]]
        end

        if @inbounds rxnarray[5,rxnind] == 0
            @fastmath @inbounds rderiv = kfs[rxnind]*krevs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[4,rxnind]]
        elseif @inbounds rxnarray[6,rxnind] == 0
            @fastmath @inbounds rderiv = kfs[rxnind]*krevs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]
        else
            @fastmath @inbounds rderiv = kfs[rxnind]*krevs[rxnind]/(p[Nspcs+rxnind]*p[Nspcs+rxnind])*cs[rxnarray[4,rxnind]]*cs[rxnarray[5,rxnind]]*cs[rxnarray[6,rxnind]]
        end

        @fastmath flux = fderiv-rderiv
        _spreadreactantpartials!(jacp,flux,rxnarray,rxnind,Nspcs+rxnind)
        _spreadproductpartials!(jacp,-flux,rxnarray,rxnind,Nspcs+rxnind)
        
        @fastmath @inbounds gderiv = rderiv*(p[Nspcs+rxnind]*p[Nspcs+rxnind])/kfs[rxnind]*RTinv

        @inbounds jacp[rxnarray[1,rxnind],rxnarray[1,rxnind]] -= gderiv
        @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[1,rxnind])
        if @inbounds rxnarray[2,rxnind] !== 0
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[1,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[1,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[2,rxnind],rxnarray[2,rxnind]] -= gderiv
            @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[2,rxnind])
            if @inbounds rxnarray[3,rxnind] !== 0
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[1,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[2,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[1,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[2,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[3,rxnind],rxnarray[3,rxnind]] -= gderiv
                @inbounds _spreadreactantpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[3,rxnind])
            end
        end

        @inbounds jacp[rxnarray[4,rxnind],rxnarray[4,rxnind]] -= gderiv
        @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[4,rxnind])
        if @inbounds rxnarray[5,rxnind] !== 0
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[4,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[4,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds jacp[rxnarray[5,rxnind],rxnarray[5,rxnind]] -= gderiv
            @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[5,rxnind])
            if @inbounds rxnarray[6,rxnind] !== 0
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[4,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[5,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[4,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[5,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds jacp[rxnarray[6,rxnind],rxnarray[6,rxnind]] -= gderiv
                @inbounds _spreadproductpartials!(jacp,gderiv,rxnarray,rxnind,rxnarray[6,rxnind])
            end
        end
    end
    jacp .*= V
end

function jacobianp!(jacp::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTPDomain,ParametrizedTPDomain}}
    jacp.=0.0
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    
    Nspcs=length(cs)
    Nrxns = size(domain.rxnarray)[2]
    
    jacobianpnsderiv!(jacp,y,p,t,domain,domain.rxnarray,cs,T,V,kfs,krevs,Nspcs,Nrxns)
    
    @simd for i in 1:length(p)
        @views @inbounds @fastmath jacp[domain.indexes[3],i] = sum(jacp[domain.indexes[1]:domain.indexes[2],i])*R*T/P
    end
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind,:] .= 0.
    end
end

function jacobianp!(jacp::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantVDomain,ParametrizedVDomain}}
    jacp.=0.0
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)

    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]
    
    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V

    jacobianpnsderiv!(jacp,y,p,t,domain,domain.rxnarray,cs,T,V,kfs,krevs,Nspcs,Nrxns)

    @simd for i in 1:Nspcs
        @views @fastmath @inbounds jacp[domain.indexes[3],i] = -(dot(Us,jacp[domain.indexes[1]:domain.indexes[2],i])+dydt[i])/(N*Cvave)
        @views @fastmath @inbounds jacp[domain.indexes[4],i] = sum(jacp[domain.indexes[1]:domain.indexes[2],i])*R*T/V + P/T*jacp[domain.indexes[3],i]
    end

    @simd for i in Nspcs+1:Nspcs+Nrxns
        @views @fastmath @inbounds jacp[domain.indexes[3],i] = -dot(Us,jacp[domain.indexes[1]:domain.indexes[2],i])/(N*Cvave)
        @views @fastmath @inbounds jacp[domain.indexes[4],i] = sum(jacp[domain.indexes[1]:domain.indexes[2],i])*R*T/V + P/T*jacp[domain.indexes[3],i]
    end

    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            @simd for i in 1:Nspcs
                ddGidTdt = flow*(-ns[i]/N)/(N*Cvave)
                jacp[domain.indexes[3],i] += ddGidTdt
                jacp[domain.indexes[4],i] += P/T*ddGidTdt
            end
        end
    end
end

function jacobianp!(jacp::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantPDomain,ParametrizedPDomain}}
    jacp.=0.0
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    
    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]
    
    dydt = zeros(size(y))
    addreactionratecontributions!(dydt,domain.rxnarray,cs,kfs,krevs)
    dydt .*= V
    
    jacobianpnsderiv!(jacp,y,p,t,domain,domain.rxnarray,cs,T,V,kfs,krevs,Nspcs,Nrxns)
    
    @fastmath Cpave = Cvave+R
    @simd for i in 1:Nspcs
        @views @fastmath @inbounds jacp[domain.indexes[3],i] = -(dot(Hs,jacp[domain.indexes[1]:domain.indexes[2],i])+dydt[i])/(N*Cpave) #divide by V to cancel ωV to ω
        @views @fastmath @inbounds jacp[domain.indexes[4],i] = sum(jacp[domain.indexes[1]:domain.indexes[2],i])*R*T/P + jacp[domain.indexes[3],i]*V/T
    end

    @simd for i in Nspcs+1:Nspcs+Nrxns
        @views @fastmath @inbounds jacp[domain.indexes[3],i] = -(dot(Hs,jacp[domain.indexes[1]:domain.indexes[2],i]))/(N*Cpave) #divide by V to cancel ωV to ω
        @views @fastmath @inbounds jacp[domain.indexes[4],i] = sum(jacp[domain.indexes[1]:domain.indexes[2],i])*R*T/P + jacp[domain.indexes[3],i]*V/T
    end
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind,:] .= 0.
    end
    
    @simd for inter in interfaces
        if isa(inter,Inlet) && domain == inter.domain
            flow = inter.F(t)
            for i in 1:Nspcs
                ddGidTdt = flow*(- ns[i]/N)/(N*Cpave)
                jacp[domain.indexes[3],i] += ddGidTdt
                jacp[domain.indexes[4],i] += ddGidTdt*V/T 
            end
        end
    end
end

function jacobianp!(jacp::Q,y::U,p::W,t::Z,domain::D,interfaces::Q3,colorvec::Q2=nothing) where {Q3<:AbstractArray,Q2,Q<:AbstractArray,U<:AbstractArray,W,Z<:Real,D<:Union{ConstantTVDomain,ParametrizedTConstantVDomain,ConstantTAPhiDomain}}
    jacp.=0.0
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(domain,y,t,p)
    
    Nspcs = length(cs)
    Nrxns = size(domain.rxnarray)[2]
    
    jacobianpnsderiv!(jacp,y,p,t,domain,domain.rxnarray,cs,T,V,kfs,krevs,Nspcs,Nrxns)
    
    @simd for ind in domain.constantspeciesinds
        @inbounds jacp[ind,:] .= 0.
    end
end

export jacobianp!

function getreactionindices(ig::Q) where {Q<:AbstractPhase}
    return deepcopy(ig.rxnarray)
end
export getreactionindices

@inline function getsensspcsrxns(domain::D,ind::Int64) where {D<:AbstractDomain}
    sensspcinds = Array{Int64,1}()
    sensrxninds = Array{Int64,1}()
    for rxnind in 1:size(domain.rxnarray)[2]
        if ind in @inbounds @view domain.rxnarray[:,rxnind]
            for spcind in @inbounds @view domain.rxnarray[:,rxnind]
                if !(spcind in sensspcinds) && (spcind !== 0)
                    push!(sensspcinds,spcind)
                end
            end
            push!(sensrxninds,rxnind)
        end
    end
    
    sensrxns = Array{ElementaryReaction,1}(undef,length(sensrxninds))
    sensspcs = Array{Species,1}(undef,length(sensspcinds))
    sensspcnames = Array{String,1}(undef,length(sensspcinds))
    senstooriginspcind = Array{Int64,1}(undef,length(sensspcinds))
    senstooriginrxnind = Array{Int64,1}(undef,length(sensrxninds))
    for (i,spcind) in enumerate(sensspcinds)
        spc = domain.phase.species[spcind]
        sensspcnames[i] = spc.name
        @inbounds sensspcs[i] = Species(
            name=spc.name,
            index=i,
            inchi=spc.inchi,
            smiles=spc.smiles,
            thermo=spc.thermo,
            atomnums=spc.atomnums,
            bondnum=spc.bondnum,
            diffusion=spc.diffusion,
            radius=spc.radius,
            radicalelectrons=spc.radicalelectrons,
        )
        @inbounds senstooriginspcind[i] = spcind
    end
    
    for (i, rxnind) in enumerate(sensrxninds)
        rxn = domain.phase.reactions[rxnind]
        reactants = Array{Species,1}()
        reactantinds = Array{Int64,1}()
        @simd for reactant in rxn.reactants
            ind = findfirst(isequal(reactant.name),sensspcnames)
            @inbounds push!(reactants,sensspcs[ind])
            push!(reactantinds,ind)
        end
        products = Array{Species,1}()
        productinds = Array{Int64,1}()
        @simd for product in rxn.products
            ind = findfirst(isequal(product.name),sensspcnames)
            @inbounds push!(products,sensspcs[ind])
            push!(productinds,ind)
        end

        @inbounds sensrxns[i] = ElementaryReaction(
            index=i,
            reactants=SVector(reactants...),
            reactantinds=MVector(reactantinds...),
            products=SVector(products...),
            productinds=MVector(productinds...),
            kinetics=rxn.kinetics,
            radicalchange=rxn.radicalchange,
            pairs=rxn.pairs
        )
        @inbounds senstooriginrxnind[i] = rxnind
    end
    
    return sensspcs,sensrxns,sensspcnames,senstooriginspcind,senstooriginrxnind
end

@inline function getsensdomain(domain::D,ind::Int64) where {D<:AbstractDomain}
    
    sensspcs,sensrxns,sensspcnames,senstooriginspcind,senstooriginrxnind = getsensspcsrxns(domain,ind)
    
    initialconds = Dict{String,Float64}()

    for fn in fieldnames(typeof(domain))
        if fn in (:T, :P, :V)
            initialconds["$fn"] = getfield(domain,fn)
        end
    end
    
    d = split(repr(typeof(domain)),"{")[1]
    
    if occursin(".",d)
        d = Symbol(split(d,".")[2])
    else 
        d = Symbol(d)
    end
    
    if isa(domain.phase,IdealGas)
        return eval(d)(phase=IdealGas(sensspcs,sensrxns,name="phase"),initialconds=initialconds)[1],sensspcnames,senstooriginspcind,senstooriginrxnind
    else
        return eval(d)(phase=IdealDiluteSolution(sensspcs,sensrxns,domain.phase.solvent,name="phase"),initialconds=initialconds)[1],sensspcnames,senstooriginspcind,senstooriginrxnind
    end
end
