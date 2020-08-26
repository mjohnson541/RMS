using Test
using DiffEqBase
using Sundials

@testset "Test Reactors" begin

phaseDict = readinput("../src/testing/liquid_phase.rms")
spcs = phaseDict["phase"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
rxns = phaseDict["phase"]["Reactions"];
solv = phaseDict["Solvents"][1];
liq = IdealDiluteSolution(spcs,rxns,solv;name="phase",diffusionlimited=true) #Define the phase (how species thermodynamic and kinetic properties calculated)

@testset "Test liquid phase Constant T Constant V reactor simulation" begin
#Constant T and V Ideal Dilute Liquid
initialconds = Dict(["T"=>450.0,"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantTVDomain(phase=liq,initialconds=initialconds,constantspecies=["oxygen"]) #Define the domain (encodes how system thermodynamic properties calculated)
react = Reactor(domain,y0,(0.0,140000.01);p=p) #Create the reactor object

sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-8); #solve the ode associated with the reactor

spcnames = getfield.(liq.species,:name)
octaneind = findfirst(isequal("octane"),spcnames)
y = sol(32977.61568)
@test y[octaneind]/sum(y) ≈ 0.461599061 rtol=3e-2 #from RMG simulator I believe the slight difference is due to better calculation of diffusion limits in RMS

#analytic jacobian vs. ForwardDiff jacobian
t=32977.61568;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j=jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;

@testset "Test liquid phase Parametrized T Constant V reactor jacobian" begin
#Parametrized T constant V Ideal Dilute Liquid
initialconds = Dict(["ts"=>[0.,600.,1200.],"T"=>[450.0,490.,500.],"V"=>1.0e-6*1e6,"octane"=>6.154e-3*1e6,"oxygen"=>4.953e-6*1e6])
domain,y0,p = ParametrizedTConstantVDomain(phase=liq,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,140000.01),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-8); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=32977.61568;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j=jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;

#Use superminimal example to test
phaseDict = readinput("../src/testing/superminimal.rms") #load mechanism dictionary
spcs = phaseDict["phase"]["Species"]; #mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
rxns = phaseDict["phase"]["Reactions"];
ig = IdealGas(spcs,rxns,name="phase")

#Constant T and P Ideal Gas
@testset "Test constant T and P reactor simulation" begin
 #Define the phase (how species thermodynamic and kinetic properties calculated)
initialconds = Dict(["T"=>1000.0,"P"=>1e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,150.11094);p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor
sim = Simulation(sol,domain);

spcnames = getfield.(ig.species,:name)
h2ind = findfirst(isequal("H2"),spcnames)
o2ind = findfirst(isequal("O2"),spcnames)
h2oind = findfirst(isequal("H2O"),spcnames)
y = sol(20.44002454)
N = sim.N(20.44002454)
@test y[h2ind]/N ≈ 0.412883111 rtol=1e-4 #from RMG simulator
@test y[o2ind]/N ≈ 0.200419093 rtol=1e-4
@test y[h2oind]/N ≈ 0.386618602 rtol=1e-4

#analytic jacobian vs. ForwardDiff jacobian
t=20.44002454;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

#sensitivities
dps = getadjointsensitivities(sim,"H2",CVODE_BDF(linear_solver=:GMRES);sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),abstol=1e-16,reltol=1e-6)
react2 = Reactor(domain,y0,(0.0,150.11094);p=p,forwardsensitivities=true)
sol2 = solve(react2.ode,CVODE_BDF(linear_solver=:GMRES),abstol=1e-21,reltol=1e-7); #solve the ode associated with the reactor
sim2 = Simulation(sol2,domain)

x,dp = extract_local_sensitivities(sol2,150.11094);
ind = findfirst(isequal("H2"),sim2.names)
dpvs = [v[ind] for v in dp]
dpvs[length(domain.phase.species)+1:end] .*= domain.p[length(domain.phase.species)+1:end]
dpvs ./= sol2(150.11094)[ind]
rerr = (dpvs .- dps')./dpvs
rerr = [isinf(x) ? 0.0 : x for x in rerr]
@test all((abs.(rerr) .> 1e-1).==false)
end;

#Constant V adiabatic Ideal Gas
#uses superminimal.yml mechanism
@testset "Constant volume adiabatic reactor simulation" begin
initialconds = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantVDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

ts = exp.(range(log(1e-15),length=10000,stop=log(0.1)))
IDT = ts[argmax(diff([sol(t)[end] for t in ts]))] #Ignition Delay Time based on argmax(dTdt(t))

@test IDT ≈ 0.038384723436228063 rtol=1e-5 #from Cantera simulation

#analytic jacobian vs. ForwardDiff jacobian
t=0.01;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;

#Constant P adiabatic Ideal Gas
#uses ethane.rms mechanism
@testset "Constant pressure adiabatic reactor simulation" begin

phaseDict = readinput("../src/testing/ethane.rms")
spcs = phaseDict["phase"]["Species"]
rxns = phaseDict["phase"]["Reactions"]
ig = IdealGas(spcs,rxns,name="phase")

initialconds = Dict(["T"=>1000.0,"P"=>2.0e5,"ethane"=>1.0,"Ar"=>1.0,"O2"=>3.5]) #Set simulation Initial Temp and Pressure
domain,y0,p = ConstantPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,1.0);p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor

ts = exp.(range(log(1e-15),length=10000,stop=log(0.2)))
IDT = ts[argmax(diff([sol(t)[end] for t in ts]))] #Ignition Delay Time based on argmax(dTdt(t))

@test IDT ≈ 0.07324954954380769 rtol=1e-5

#analytic jacobian vs. ForwardDiff jacobian
t=0.01;
y=sol(t)
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)

end;

#Parametrized T and P Ideal Gas
#uses superminimal.yml mechanism
@testset "Parametrized T and P reactor jacobian" begin
initialconds = Dict(["T"=>[1000.0,1400.0,1500.0],"ts"=>[0.,100.,200.],"P"=>[1.0e5,1.8e5,2.0e5],"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp/Pressure and Volume (function/array)
domain,y0,p = ParametrizedTPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;

#Parametrized V adiabatic Ideal Gas
#uses superminimal.yml mechanism
@testset "Parametrized volume adiabatic reactor jacobian" begin
initialconds = Dict(["T"=>1000.0,"ts"=>[0.,100.,200.],"V"=>[0.01,0.05,0.1],"H2"=>0.67,"O2"=>0.33]) #Set simulation Initial Temp/Pressure and Volume (function/array)
domain,y0,p = ParametrizedVDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;

#Parametrized P adiabatic Ideal Gas
#uses ethane.yml mechanism
@testset "Parametrized pressure adiabatic reactor jacobian" begin
initialconds = Dict(["T"=>1000.0,"ts"=>[0.,100.,200.],"P"=>[2.0e5,3.0e5,5.0e5],"H2"=>0.67,"O2"=>0.33])
domain,y0,p = ParametrizedPDomain(phase=ig,initialconds=initialconds) #Define the domain (encodes how system thermodynamic properties calculated)

react = Reactor(domain,y0,(0.0,0.101),p=p) #Create the reactor object
sol = solve(react.ode,CVODE_BDF(),abstol=1e-20,reltol=1e-12); #solve the ode associated with the reactor

#analytic jacobian vs. ForwardDiff jacobian
t=0.1;
y = sol(t);
ja=jacobiany(y,p,t,domain,[],nothing);
j = jacobianyforwarddiff(y,p,t,domain,[],nothing);
@test all((abs.(ja.-j) .> 1e-4.*abs.(j).+1e-16).==false)
end;



@testset "Multi-domain ConstantV sensitivity analysis" begin
    phaseDict = readinput("../src/testing/superminimal.rms")
    spcs = phaseDict["phase"]["Species"]
    rxns = phaseDict["phase"]["Reactions"]
    ig = IdealGas(spcs,rxns,name="phase")
    
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactV = Reactor(domainV,y0V,(0.0,0.037);p=pV) #Create the reactor object
    solV = solve(reactV.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    simV = Simulation(solV,domainV)
    
    initialcondsV1 = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV1,y0V1,pV1 = ConstantVDomain(phase=ig,initialconds=initialcondsV1) #Define the domain (encodes how system thermodynamic properties calculated)
    initialcondsV2 = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV2,y0V2,pV2 = ConstantVDomain(phase=ig,initialconds=initialcondsV2) #Define the domain (encodes how system thermodynamic properties calculated)
    
    react,y0,p = Reactor((domainV1,domainV2),(y0V1,y0V2),(0.0,0.037),[],(pV1,pV2));
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6);
    sysim = SystemSimulation(sol,(domainV1,domainV2),p);
    
    t = 0.03
    @test sol(t)[1:length(spcs)] ≈ solV(t)[1:end-2] rtol=1e-5
    @test sol(t)[length(spcs)+1:end-4] ≈ solV(t)[1:end-2] rtol=1e-5
    
    dpsV = getadjointsensitivities(simV,"H2",CVODE_BDF();sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),abstol=1e-16,reltol=1e-6)
    dps = getadjointsensitivities(sysim,sysim.sims[1],"H2",CVODE_BDF();sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)),abstol=1e-16,reltol=1e-6)
    @test dpsV ≈ dps rtol=1e-4
end;

@testset "Multi-domain ConstantV and ConstantTP simulation" begin
    phaseDict = readinput("../src/testing/superminimal.rms")
    spcs = phaseDict["phase"]["Species"]
    rxns = phaseDict["phase"]["Reactions"]
    ig = IdealGas(spcs,rxns,name="phase")
    
    initialcondsTP = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainTP,y0TP,pTP = ConstantTPDomain(phase=ig,initialconds=initialcondsTP) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactTP = Reactor(domainTP,y0TP,(0.0,0.04);p=pTP) #Create the reactor object
    solTP = solve(reactTP.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    reactV = Reactor(domainV,y0V,(0.0,0.04);p=pV) #Create the reactor object
    solV = solve(reactV.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6); #solve the ode associated with the reactor
    
    initialcondsTP = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainTP,y0TP,pTP = ConstantTPDomain(phase=ig,initialconds=initialcondsTP) #Define the domain (encodes how system thermodynamic properties calculated)
    initialcondsV = Dict(["T"=>1000.0,"P"=>10.0e5,"H2"=>0.67,"O2"=>0.33]) 
    domainV,y0V,pV = ConstantVDomain(phase=ig,initialconds=initialcondsV) #Define the domain (encodes how system thermodynamic properties calculated)
    
    react,y0,p = Reactor((domainTP,domainV),(y0TP,y0V),(0.0,0.04),[],(pTP,pV));
    sol = solve(react.ode,CVODE_BDF(),abstol=1e-16,reltol=1e-6);
    
    t = 0.03
    @test sol(t)[1:length(spcs)] ≈ solTP(t)[1:end-1] rtol=1e-5
    print(length(sol(t)[length(spcs)+1:end-3]))
    print(length(solV(t)[1:end-2]))
    @test sol(t)[length(spcs)+1:end-3] ≈ solV(t)[1:end-2] rtol=1e-5
end;
end;