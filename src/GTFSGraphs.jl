module GTFSGraphs

using CSV, DataFrames
using NearestNeighbors
import Distances: Metric
using Graphs, MetaGraphsNext

#include("../../GTFSQuerys/src/GTFSQuerys.jl")
using GTFSQuerys

export GTFSGraph, AbstractVertexData, AbstractEdgeData, CoordTuple, XYTuple
export add!, add_headway_and_dwell!, add_running_time!, add_s2s_time!
export apply_vertex_function!, apply_edge_function!, connect_vertexs!
export make_tables, save_graph


"""
    Abstract struct of data to store at each vertex
    Concrete subtypes must be mutable and need to implament:
    - `VertexData()` Default constructor
    - All fields of L LabelType
    - `dwell`<:Real
    - `headway`<:Real
    
"""
abstract type AbstractVertexData end

"""
    Abstract struct of data to store for each edge

    Concrete subtypes must be mutable and need to implement:
    - `EdgeData()` Default constructor
    - `running_time`<:Real time delta between stops in the schedule
    - `s2s_time`<:Real
    - `shape_distance`<:AbstractFloat distance of edge along shape
    - `edge_type`::String
"""
abstract type AbstractEdgeData end


const CoordTuple = NamedTuple{(:lat, :lon), NTuple{2, Float64}}
const XYTuple = NamedTuple{(:x, :y), NTuple{2, Float64}}
# Base.string(t::CoordTuple) = "($(t[:lat]),$(t[:lon]))"
# Base.string(t::XYTuple) = "($(t[:s]),$(t[:y]))"


mutable struct GTFSGraph{L,V<:AbstractVertexData,E<:AbstractEdgeData}
    g::MetaGraph
    q::GTFSQuery

    # Paramter sets
    shape_params::Vector{Symbol}
    vertex_label_params::Vector{Symbol}
    
    # Funciton definitions
    distance_function::Metric
end

#
# Constructors
#

function GTFSGraph{L,V,E}(query::GTFSQuery; 
                          shape_params::Vector{Symbol}=[:shape_id],
                          #vertex_label_params::Vector{Symbol}=[:stop_id, :route_id, :direction_id],
                          #vertex_data_params::Vector{Symbol}=[:stop_id, :route_id, :direction_id, :stop_lat, :stop_lon],
                          distance_function::Metric=Haversine(6371000)
                         ) where {L, V<:AbstractVertexData, E<:AbstractEdgeData}

    # TODO Checks
    # vertex_label_params == propertynams(L)
    # Feildnames of V contains vertex_label_params and :stop_coord
    # Feildnemas of E contains [:weights, :shape_distance]
    # default constructors exist for L, V, E


    vertex_label_params = collect(fieldnames(L))

    g = MetaGraph(DiGraph(),L,V,E,"Graph of Transit Network",e_data::E -> e_data.weight,Inf)

    GTFSGraph{L,V,E}(g,query,shape_params,vertex_label_params,distance_function)
end

#
# Methods
#

function add!(gg::GTFSGraph, feature_filter::FilterDict)

    # Get unique set of the shape_params
    #filters = FilterDict(:route_id => x::String -> ==(x,route_id))
    shape_df = query(gg.q, [gg.shape_params..., :trip_headsign], feature_filter)  # should apply date filters, ideally not required

    # TODO fix the unique check
    # if shape_df == unique(shape_df)
    #     throw(ArgumentError("The shape_params of $(gg.shape_params) resulted in non-unique shapes to be added to the graph"))
    # end
    
    # Loop over shape params
    for row in eachrow(shape_df)
        # TODO switch to using a grouped dataframe to minimize query calls

        # Create Filters
        filters = FilterDict()
        for p in gg.shape_params
            merge!(filters, FilterDict(p => x -> .==(x,row[p])))
        end

        # Get table of data for vertex creation
        vertex_params = unique(vcat(gg.vertex_label_params, [:stop_lat, :stop_lon, :stop_sequence]))
        v_df = query(gg.q, vertex_params, filters)  # Do not need the date filters on this one
        sort!(v_df, [order(:stop_sequence, by=x->parse(Int64,x))])

        # Get table of data for edge creation
        shape_params = [:shape_pt_lat, :shape_pt_lon, :shape_pt_sequence]
        s_df = query(gg.q, shape_params, filters)  # do not need the date filters
        sort!(s_df, [order(:shape_pt_sequence, by=x->parse(Int64,x))])

        #println(row)
        #println(v_df[:,[:stop_id, :route_id, :direction_id, :stop_name, :stop_sequence, :trip_headsign]])
        add_shape!(gg, v_df, s_df)
    end
    return nothing
end

function add_shape!(gg::GTFSGraph{L, V, E}, v_df::DataFrame, s_df::DataFrame) where {L, V, E}
    # Compute nearest shape point 
    shape_x = parse.(Float64, s_df[:,:shape_pt_lat])
    shape_y = parse.(Float64, s_df[:,:shape_pt_lon])
    tree = BallTree(vcat(shape_x', shape_y'), gg.distance_function)

    stop_x = parse.(Float64, v_df[:,:stop_lat])
    stop_y = parse.(Float64, v_df[:,:stop_lon])
    near_idx, dist = nn(tree, vcat(stop_x',stop_y'))

    # Loop over rows of v_df
        #add_vertex
    prv_label = L(["" for _ in fieldnames(L)])
    for (i, row) in enumerate(eachrow(v_df))
        #
        # Add Vertex
        #

        # Create label tuple and add empty vertex if key does not exist
        label = L(row[gg.vertex_label_params])
        !haskey(gg.g, label) ? gg.g[label] = V() : nothing    

        # Add all label data to vertex data
        for (sym, val) in zip(gg.vertex_label_params, row[gg.vertex_label_params])
            setfield!(gg.g[label], sym, val) 
        end
        setfield!(gg.g[label], :stop_coord, CoordTuple((stop_x[i],stop_y[i])))

        # Continue if on first iteration
        i == 1 ? (prv_label = label ; continue) : nothing


        #
        # Add Edge
        #
        # Compute distance between points in shape
        dist = 0.0
        for idx in near_idx[i-1]:near_idx[i]-1
            a = vcat(shape_x[idx],     shape_y[idx])
            b = vcat(shape_x[idx + 1], shape_y[idx + 1])
            dist += gg.distance_function(a,b)
        end

        # TODO add distance check to stop maybe above

        # Create edge if one does not exist
        !haskey(gg.g, prv_label, label) ? gg.g[prv_label, label] = E() : nothing

        # Set fields
        setfield!(gg.g[prv_label, label], :shape_distance, dist)
        setfield!(gg.g[prv_label, label], :edge_type, "route_segment")
        
        prv_label = label
    end
    return nothing
end

#
# Includes
#

include("utilities.jl")
include("apply_functions.jl")
include("methods.jl")

end # module GTFSGraphs

# function add_route!(gg::GTFSGraph, route_id::String)

#     # Get unique set of the shape_params
#     filters = FilterDict(:route_id => x -> .==(x,route_id))
#     shape_df = query(gg.q, gg.shape_params, filters)  # should apply date filters, ideally not required

#     # TODO add check that they are actually unique
    
#     # Loop over shape params
#     for row in eachrow(shape_df)

#         # Create Filters
#         filters = FilterDict()
#         for p in gg.shape_params
#             merge!(filters, FilterDict(p => x -> .==(x,row[p])))
#         end

#         # Get table of data for vertex creation
#         vertex_params = unique(vcat(gg.vertex_label_params, gg.vertex_data_params, [:stop_lat, :stop_lon, :stop_sequence]))
#         v_df = query(gg.q, vertex_params, filters)  # Do not need the date filters on this one
#         sort!(v_df, [order(:stop_sequence, by=x->parse(Int64,x))])

#         # Get table of data for edge creation
#         shape_params = [:shape_pt_lat, :shape_pt_lon, :shape_pt_sequence]
#         s_df = query(gg.q, shape_params, filters)  # do not need the date filters
#         sort!(s_df, [order(:shape_pt_sequence, by=x->parse(Int64,x))])

#         #println(row)
#         #println(v_df[:,[:stop_id, :route_id, :direction_id, :stop_name, :stop_sequence, :trip_headsign]])
#         add_shape!(gg, v_df, s_df)
#     end
#     return nothing
# end


# function query_stop_times(gg::GTFSGraph{L, V, E}, valid_labels::Vector{L}) where {L, V, E}

#     # Create and fill set tuple
#     prop_names = fieldnames(L)
#     set_dict = Dict{Symbol,Set{String}}(prop_names .=> [Set{String}() for _ in prop_names])
#     #set_tuple = NamedTuple{prop_names, NTuple{length(prop_names),Set{String}}}([Set{String}() for _ in prop_names])

#     for label in valid_labels
#         for prop in prop_names
#             push!(set_dict[prop],label[prop])
#         end
#     end

#     # Create Filters
#     params = [:arrival_time, :departure_time, gg.vertex_label_params...]
#     filters = FilterDict()
#     for p in gg.vertex_label_params
#         merge!(filters, FilterDict(p => x -> in.(x, [set_dict[p]])))
#     end
#     df = query(gg.q, params, filters)
#     gdf = groupby(df, gg.vertex_label_params)
#     return gdf
# end
