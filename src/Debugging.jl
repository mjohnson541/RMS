using Statistics

"""
This function prints a report analyzing rates and dydt to determine reactions
and thermochemistry worth looking at
sim should be a simulation object created from a solution object corresponding to a crash
tol is the ratio relative to the median absolute rate or dn_i/dt value above which reactions
and species will be reported
"""
function analyzecrash(sim::Simulation;tol=1e6)
    t = sim.sol.t[end]
    rts = rates(sim,t)
    rmedian = median(abs.([rt for rt in rts if !isnan(rt) && rt != 0.0]))
    dydt = zeros(length(sim.sol.u[end]))
    dydtreactor!(dydt,sol.u[end],0.0,sim.domain,[];p=sim.domain.p)
    dymedian = median(abs.([dy for dy in dydt if !isnan(dy) && dy != 0.0]))
    y = sim.sol(t)
    p = sim.domain.p
    ns,cs,T,P,V,C,N,mu,kfs,krevs,Hs,Us,Gs,diffs,Cvave,cpdivR = calcthermo(sim.domain,y,t,p)
    println("Crash Analysis Report")
    for (i,rt) in enumerate(rts)
        if isnan(rt)
            println("NaN for Reaction:")
            ind = findfirst(x->x==rxn,sim.domain.phase.reactions)
            analyzereaction(sim.domain.phase.reactions[i],sim,t,kfs[ind],krevs[ind])
        elseif abs(rt/rmedian) > tol
            ratio = abs(rt/rmedian)
            println("rt/rmedian > $tol, rt=$rt, rt/rmedian=$ratio \nReaction:")
            ind = findfirst(x->x==rxn,sim.domain.phase.reactions)
            analyzereaction(sim.domain.phase.reactions[i],sim,t,kfs[ind],krevs[ind])
        end
    end
    for (i,dy) in enumerate(dydt)
        if i > length(sim.domain.phase.species)
            continue
        end
        if isnan(dy)
            println("NaN for Species net flux:")
            analyzespecies(sim.domain.phase.species[i],sim,t)
        elseif abs(dy/dymedian) > tol
            ratio = abs(dy/dymedian)
            println("dydt/dydtmedian > $tol, dydt=$dy, dydt/dydtmedian=$ratio \nSpecies:")
            analyzespecies(sim.domain.phase.species[i],sim,t)
        end
    end
end

function analyzereaction(rxn,sim,t,kf,krev)
    println(getrxnstr(rxn))
    if length(rxn.reactants) == 1
        kunitsf = "s^-1"
    elseif length(rxn.reactants) == 2
        kunitsf = "m^3/(mol*s)"
    elseif length(rxn.reactants) == 3
        kunitsf = "m^6/(mol^2*s)"
    else
        error("cannot accomodate reactants>3")
    end
    if length(rxn.products) == 1
        kunitsr = "s^-1"
    elseif length(rxn.products) == 2
        kunitsr = "m^3/(mol*s)"
    elseif length(rxn.products) == 3
        kunitsr = "m^6/(mol^2*s)"
    else
        error("cannot accomodate products>3")
    end
    println("kf: $kf $kunitsf")
    println("krev: $krev $kunitsr")
    println("Reactants: ")
    for reactant in rxn.reactants
        analyzespecies(reactant,sim,t)
    end
    println("Products: ")
    for product in rxn.products
        analyzespecies(product,sim,t)
    end
end

function analyzespecies(spc,sim,t)
    H298 = getEnthalpy(spc.thermo,298.0)/4184.0
    name = spc.name
    println("\t$name H298: $H298 kcal/mol")
end

"""
This calculates the collision limit of H + H -> H2 as an upper bound
"""
function calccollisionlimit(T)
    mu = 0.00050397
    sigma = 2.05e-10
    epsilon = 1205.6
    Tr = T*kB*Na/epsilon
    collintegral = 1.16145*Tr^(-0.14874)+0.52487*exp(-0.7732 * Tr)+2.16178*exp(-2.437887*Tr)
    kcoll = sqrt(8.0*pi*kB*T*Na/mu)*sigma^2*collintegral*Na
    return kcoll
end

"""
Compares bimolecular and trimolecular reactions with the collision limit for H+H->H2 
to check whether they are physical, a report with only a title and kcoll means no violators
"""
function analyzecolllimit(phase,Tmin,Tmax,Pmin,Pmax)
    println("Collision Limit Report (Comparisons are done with collision limit for H+H->H2)")
    println("kcoll=$kcoll")
    for T in [Tmin,Tmax]
        for P in [Pmin,Pmax]
            kcoll = calccollisionlimit(T)
            cpdivR,hdivRT,sdivR = calcHSCpdless(phase.vecthermo,T)
            Gs = (hdivRT.-sdivR)*(R*T)
            if sim.domain.phase.diffusionlimited
                diffs = getfield.(sim.domain.phase.species,:diffusion)(T=T,mu=0.0,P=P)
            else
                diffs = Array{Float64,1}()
            end
            kfs,krevs = getkfkrevs(phase,T,P,P/(R*T),1.0,Gs,diffs,V=1.0/C)
            for (i,rxn) in enumerate(phase.reactions)
                boo = false
                if len(rxn.reactants) > 2
                    if kfs[i] > kcoll
                        if !boo
                            println(getrxnstr(rxn))
                            boo = true
                        end
                        kf = kfs[i]
                        ratio = kfs[i]/kcoll
                        println("Violated in forward direction kf=$kf ratio=$ratio T=$T P=$P")
                    end
                end
                if len(rxn.products) > 2
                    if krevs[i] > kcoll
                        if !boo
                            println(getrxnstr(rxn))
                            boo = true
                        end
                        krev = krevs[i]
                        ratio = krevs[i]/kcoll
                        println("Violated in forward direction krev=$krev ratio=$ratio T=$T P=$P")
                    end
                end
            end
        end
    end
end