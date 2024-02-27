# function returns grid with required face-/nodesets
function load_model(modeltype::String)
    path = string("mesh/", modeltype,".msh")
    grid = FerriteGmsh.togrid(path)
    
    # boundary box
    addfaceset!(grid, "left", x -> x[1] ≈ 0.0)
    addfaceset!(grid, "right", x -> x[1] ≈ 0.09)
    addfaceset!(grid, "top", x -> x[2] ≈ 0.0)
    addfaceset!(grid, "bottom", x -> x[2] ≈ 0.09)
    addfaceset!(grid, "front", x -> x[3] ≈ 0.0)
    addfaceset!(grid, "back", x -> x[3] ≈ 0.36)
    
    # cylinder surfaces
    addfaceset!(grid, "leftsurf", x -> x[3] ≈ 0.03 && sqrt((x[1]-0.045)^2+(x[2]-0.045)^2) <= 0.016) 
    addfaceset!(grid, "rightsurf", x -> x[3] ≈ 0.33 && sqrt((x[1]-0.045)^2+(x[2]-0.045)^2) <= 0.016)
    addfaceset!(grid, "middlesurf", difference(in_rectangle([0.03,0.06,0.18],[0.04,0.04,0.05]),in_cylinder([0.045,0.045,0.18],[0.0,0.0,1.0],0.013,0.3)))
    
    # single nodes
    addnodeset!(grid, "luv", x -> x[1] ≈ 0.0 && x[2] ≈ 0.0 && x[3] ≈ 0.0)
    addnodeset!(grid, "luh", x -> x[1] ≈ 0.0 && x[2] ≈ 0.0 && x[3] ≈ 0.36)
    addnodeset!(grid, "rov", x -> x[1] ≈ 0.09 && x[2] ≈ 0.09 && x[3] ≈ 0.0)
    
    # surrounding medium
    addfaceset!(grid, "airbc", difference(in_rectangle([0.045,0.045,0.18],[0.1,0.1,0.37]),in_cylinder([0.045,0.045,0.18],[0.0,0.0,1.0],0.016,0.305)))     
    
    return grid
end