"""
Returns only free surfaces of given faceset
"""
function getsurfaceset(grid::Grid,faces::Set{Tuple{Int64,Int64}})
    surfset = Set{Tuple{Int,Int}}()
    for face in faces
        if grid.boundary_matrix[face[2],face[1]][1]
            push!(surfset,face)
        end
    end
    return surfset
end

"""
Routines return functions to check if a given point lies within specific domains defined by basic geometries.

Routines also except 0.0 lengths and Inf.

Tolerance tol can be set to negative value to exclude points on domain border

Multiple geometry sets can be combined with 'intersecting' (logical &&) and 'unifying' (logical ||). Example:
 'unifying(in_point([0.0]),in_point([1.0]))' -> true for x ∈ {[0.0],[1.0]}
 'intersecting(beyond_basis_plane(1,1.0),in_sphere([1.0,0.0],1.0))'' -> true for x ∈ right half of cylinder with radius 1.0 and mid-point [1.0,0.0]
"""

function intersecting(f1::Function,f2::Function)
    return x -> f1(x) && f2(x)
end

function unifying(f1::Function,f2::Function)
    return x -> f1(x) || f2(x)
end

function difference(f1::Function,f2::Function)
    return x -> f1(x) && ! f2(x)
end

function in_point(point::Vector{Float64};tol::Float64=1.0e-6)
    return x -> norm(x-point) ≤ abs(tol)
end


function in_sphere(mid::Vector{Float64},radius::Float64;tol::Float64=1.0e-6)
    if radius == 0.0 # actually point
        return in_point(mid;tol=tol)
    else
        return x -> norm(x-mid) ≤ radius+tol
    end
end

function beyond_plane(basePoint::Vector{Float64},normal::Vector{Float64};tol::Float64=1.0e-6)
    return x -> basePoint ⋅ normal - x ⋅ normal ≥ tol
end

"""
apply negativ basis direction with dir < 0
"""
function beyond_basis_plane(dir::Int,dist::Float64;tol::Float64=1.0e-6)
    return x -> sign(dir)*(x[abs(dir)] - dist) ≥ tol
end

function in_plane(basePoint::Vector{Float64},normal::Vector{Float64};tol::Float64=1.0e-6)
    return x -> abs(basePoint ⋅ normal - x ⋅ normal) ≤ abs(tol)
end

function in_basis_plane(dir::Int,dist::Float64;tol::Float64=1.0e-6)
    return x -> abs(x[dir] - dist) ≤ abs(tol)
end

function in_rectangle(midPoint::Vector{Float64},normal::Vector{Float64},dimensions::Vector{Float64};tol::Float64=1.0e-6)
    # 1 Normal -> 2D
    otho = [-normal[2],normal[1]]
    return x -> (
        (x - midPoint) ⋅ normal/norm(normal)  ≥ - 0.5*dimensions[1]-tol
    &&  (x - midPoint) ⋅ normal/norm(normal)  ≤   0.5*dimensions[1]+tol
    &&  (x - midPoint) ⋅ otho/norm(otho)  ≥ - 0.5*dimensions[2]-tol
    &&  (x - midPoint) ⋅ otho/norm(otho)  ≤   0.5*dimensions[2]+tol
    )
end

function in_rectangle(midPoint::Vector{Float64},normal1::Vector{Float64},normal2::Vector{Float64},dimensions::Vector{Float64};tol::Float64=1.0e-6)
    # 2 Normal -> 3D
    normal3 = normal1 × normal2
    return x -> (
        (x - midPoint) ⋅ normal1/norm(normal1)  ≥ - 0.5*dimensions[1]-tol
    &&  (x - midPoint) ⋅ normal1/norm(normal1)  ≤   0.5*dimensions[1]+tol
    &&  (x - midPoint) ⋅ normal2/norm(normal2)  ≥ - 0.5*dimensions[2]-tol
    &&  (x - midPoint) ⋅ normal2/norm(normal2)  ≤   0.5*dimensions[2]+tol
    &&  (x - midPoint) ⋅ normal3/norm(normal3)  ≥ - 0.5*dimensions[3]-tol
    &&  (x - midPoint) ⋅ normal3/norm(normal3)  ≤   0.5*dimensions[3]+tol
    )
end

function in_rectangle(midPoint::Vector{Float64},dimensions::Vector{Float64};tol::Float64=1.0e-6)
    if size(dimensions,1) == 2
        return in_rectangle(midPoint,[1.0,0.0],dimensions,tol)
    elseif size(dimensions,1) == 3
        return in_rectangle(midPoint,[1.0,0.0,0.0],[0.0,1.0,0.0],dimensions,tol=tol)
    else
        error("Rectangle dimensions are not compatible!")
        return nothing
    end
end

function in_cylinder(midPoint::Vector{Float64},normal::Vector{Float64},
                    radius::Float64,height::Float64;tol::Float64=1.0e-6)
    if size(midPoint,1) ≠ 3 || size(normal,1) ≠ 3
        error("Function 'in_cylinder()' only works with 3D-basis!")
    end
    if radius == Inf && height == 0.0 # actually plane
        return in_plane(midPoint,normal,tol=tol)
    elseif radius == 0.0 && height == 0.0 # actually point
        return in_point(midPoint,tol=tol)
    else # cylinder
        return x -> (
                (x - midPoint) ⋅ normal/norm(normal)  ≥ - 0.5*height-tol  # "above bottom"
        &&      (x - midPoint) ⋅ normal/norm(normal)  ≤   0.5*height+tol  # "below top"
        && norm((midPoint - x) × normal/norm(normal)) ≤   radius+tol      # within radius
        )  # within radius
    end
end