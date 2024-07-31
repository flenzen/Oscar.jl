# Generates the coboundary matrices of a diameter-filtered Vietoris-Rips complex.
# A Vietoris-Rips complex on a point cloud X is the abstract simplicial complex
# VR(X) ≔ 2^X. It is filtered by the diameter diam σ ≔ max{d(x, y) | x, y ∈ σ}.
# In practise, one truncates VR(X) at a fixed threshold t, i.e.,
# VR(X) ≔ {σ ∈ 2^X | diam σ ≤ t}.

#TODO Use UInt instead of Int, where applicable.
#TODO Take diameter into consideration

include("CombinatorialNumberSystem.jl")

""" Represents one stage during the computation of Vietoris-Rips coboundary maps."""
struct VietorisRipsComplex
	distances::Array{Float64, 2}
	neighbors::Vector{Vector{Int}}
	dimension::Int
	simplices::Vector{Int}
	diameters::Vector{Int}
end

"""
	vietoris_rips_complex(distances[, threshold])
	
Creates a new Vietoris-Rips complex. """
function vietoris_rips_complex(distances::Array{Float64, 2}, threshold = Inf::Float64)
	@req is_symmetric(distances) && all(distances .>= 0) "Matrix must be a distance matrix"
	n = size(distances, 1)
	VietorisRipsComplex(
		distances,
		[[j for j in 1:n if distances[i, j] <= threshold && i != j] for i in 1:n],
		# [sort([j for j in 1:n if distances[i, j] <= threshold && i != j], by=j->distances[i,j]) for i in 1:n],
		0,
		collect(1:n),
		[0 for _ in 0:n-1],
	)
end


"""
	coboundary_map(K :: VietorisRipsComplex)
	
Computes the next coboundary map of a Vietoris-Rips complex.

Returns the next piece of data, plus the actual coboundary matrix.

# Examples
```jldoctest
julia> distance_matrix = [0.0 1.0 2.0 3.0; 1.0 0.0 1.0 2.0; 2.0 1.0 0.0 1.0; 3.0 2.0 1.0 0.0];

julia> V = vietoris_rips_complex(distance_matrix);

julia> V, δ0 = coboundary_map(V); δ0
[ 1    1    0    1    0    0]
[-1    0    1    0    1    0]
[ 0   -1   -1    0    0    1]
[ 0    0    0   -1   -1   -1]

julia> V, δ1 = coboundary_map(V); δ1
[ 1    1    0    0]
[-1    0    1    0]
[ 1    0    0    1]
[ 0   -1   -1    0]
[ 0    1    0   -1]
[ 0    0    1    1]
```"""
function coboundary_map(K::VietorisRipsComplex)
	# Algorithm taken from https://gitlab.com/flenzen/2pac/-/blob/main/CliqueComplex.cpp
	coboundary_columns = [Vector{Tuple{Float64, Int, Int}}() for _ in 1:length(K.simplices)]
	simplices_next_dim = Vector{Int}()
	diameters_next_dim = Vector{Float64}()
	index_of_simplex_next_dim = Dict{Int, Int}()

	for (column_index, (σ, diam)) in enumerate(zip(K.simplices, K.diameters))
		if K.dimension == 0
			# σ = {v} = v is a vertex.
			for n in K.neighbors[σ]
				coface = encode(sort([σ, n]))
				diameter = K.distances[σ, n]
				push!(coboundary_columns[column_index], (diameter, coface, n < σ ? -1 : 1))
			end
		else
			# σ = encode(xs) is a higher dimensional simplex
			xs = decode(σ, K.dimension)
			# Iterate through neighbors of ith and i+1st vertex in xs in parallel until a common neighbour is found.
			# If indices[i] has been incremented, it is out of sync with neighbor of i-1st vertex. Therefore, continue with i-1 and i.
			# Otherwise, continue with i+1, i+2. Eventually, we have a common neighbor of all vertices in σ,
			# viz the apex of a coface.
			neighbor_indices = [1 for _ in 1:length(xs)]
			i = 1
			while neighbor_indices[i] <= length(K.neighbors[xs[i]]) && neighbor_indices[i+1] <= length(K.neighbors[i+1])
				if K.neighbors[xs[i]][neighbor_indices[i]] < K.neighbors[xs[i+1]][neighbor_indices[i+1]]
					# The i-th neighbor index is behind
					neighbor_indices[i] += 1
					if i > 1
						i -= 1
					end
				elseif K.neighbors[xs[i]][neighbor_indices[i]] > K.neighbors[xs[i+1]][neighbor_indices[i+1]]
					neighbor_indices[i+1] += 1
				elseif (i += 1) == length(xs)
					# All indices[1],...,indices[i] point to the same common neighbor of all vertices xs,
					# so we have found a coface of σ
					apex = K.neighbors[xs[i]][neighbor_indices[i]]
					diameter = maximum(K.distances[v, apex] for v in xs; init = diam)
					coface = encode(sort([xs; [apex]]))
					coeff = 1
					for v in xs
						if v > apex
							break
						else
							coeff *= -1
						end
					end
					push!(coboundary_columns[column_index], (diameter, coface, coeff))
					# Continue with finding more cofaces
					i = 1
					neighbor_indices[i] += 1
				end
			end
			# At this point, there are no more cofaces
		end
	end

	new_simplices = sort(collect(Set(σ for column in coboundary_columns for (_, σ, _) in column)))

	#TODO Replace by (heap)-sparse matrix
	dense_mat = zero_matrix(QQ, length(K.simplices), length(new_simplices))
	for (j, column) in enumerate(coboundary_columns)
		for (_, i, coeff) in column
			dense_mat[j, i+1] = coeff
		end
	end

	return VietorisRipsComplex(
		K.distances,
		K.neighbors,
		K.dimension + 1,
		new_simplices,
		[0.0 for _ in 1:length(new_simplices)],
	), dense_mat
end

