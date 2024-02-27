# applies macroscopic boundary conditions to mesh
function create_bc_macro(dh::DofHandler, amplitude::Float64)
    dbc = ConstraintHandler(dh)
    
    # fix cylinder bases
    add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "leftsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))
    add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "rightsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))

    # add mechanical displacement to middle of cylinder as main deformation driving force
    add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "middlesurf"), (x,t) -> amplitude*time_magnitude(Float64(t), n_timesteps), [1]))
    
    # add electric / magnetic potential grounding
    add!(dbc, Dirichlet(:phi, getnodeset(dh.grid, "luv"), (x,t) -> zero(Vec{1}), [1]))
    
    add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "luv"), (x,t) -> zero(Vec{3}), [1,2,3]))
    add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "luh"), (x,t) -> zero(Vec{3}), [1,2,3]))
    add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "rov"), (x,t) -> zero(Vec{3}), [1,2,3]))

    # set mechanical displacement to zero for surrounding medium air
    add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "airbc"),  (x,t) -> zero(Vec{3}), [1,2,3]))

    close!(dbc)
    t = 0.0 # t = time
    Ferrite.update!(dbc, t)
    return dbc
end 

# returns amplitude of mechanical displacement depending on timestep
function time_magnitude(time::Float64, maxtime::Int) # 0 -> 1 -> -1 -> 0 -cycle - maxtime-1 needs to be perfectly divisible by 4
    @assert mod(maxtime-1,4)==0
    t = Int(time)
    al = Int((maxtime-1)/4+1)
    bl = Int((maxtime-1)/2+1)
    cl = Int((maxtime-1)/4+1)
    a = range(0.0,1.0,length=al)
    b = range(1.0,-1.0,length=bl)
    c = range(-1.0,0.0,length=cl)
    if(t==0)
        f_m = 0.0
    elseif(t < al)
        f_m = a[t]
    elseif(t < (al+bl-1))
        f_m = b[t-al+1]
    else
        f_m = c[t-(al+bl-2)]
    end
    return f_m 
end